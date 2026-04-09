#pragma once

#include <volt/analysis/cna_classifier.h>
#include <volt/analysis/nearest_neighbor_finder.h>
#include <volt/topology/crystal_coordination_pattern.h>

#include <array>
#include <cmath>

namespace Volt::CnaLocalStructureUtils{

enum class NativeCnaLocalEnvironmentConstruction{
    None = 0,
    Direct,
    FirstShellExpansion
};

struct NativeCnaLocalEnvironmentDescriptor{
    NativeCnaLocalEnvironmentConstruction construction = NativeCnaLocalEnvironmentConstruction::None;
    int expansionSeedCount = 0;
    int referenceNeighborOffset = 0;
    int referenceNeighborCount = 0;
    double cutoffMultiplier = 0.0;
    int bondStartIndex = 0;
    int extraNeighborRejectIndex = -1;
};

struct LocalStructureMatch{
    double localCutoff = 0.0;
    int coordinationNumber = 0;
    CoordinationStructureType coordinationType = COORD_OTHER;
    StructureType structureType = StructureType::OTHER;
    std::array<int, MAX_NEIGHBORS> neighborIndices{};
    std::array<int, MAX_NEIGHBORS> neighborMapping{};
    std::array<Vector3, MAX_NEIGHBORS> neighborVectors{};
};

inline NativeCnaLocalEnvironmentDescriptor nativeCnaLocalEnvironmentFor(
    LatticeStructureType inputCrystalType
){
    switch(inputCrystalType){
        case LATTICE_FCC:
        case LATTICE_HCP:
            return {
                NativeCnaLocalEnvironmentConstruction::Direct,
                0,
                0,
                12,
                1.2071067811865475,
                0,
                12
            };
        case LATTICE_BCC:
            return {
                NativeCnaLocalEnvironmentConstruction::Direct,
                0,
                0,
                8,
                1.393846850117352,
                0,
                14
            };
        case LATTICE_CUBIC_DIAMOND:
        case LATTICE_HEX_DIAMOND:
            return {
                NativeCnaLocalEnvironmentConstruction::FirstShellExpansion,
                4,
                4,
                12,
                1.2071068,
                4,
                -1
            };
        default:
            return {};
    }
}

inline int coordinationNumberFor(LatticeStructureType inputCrystalType){
    switch(inputCrystalType){
        case LATTICE_FCC:
        case LATTICE_HCP:
            return 12;
        case LATTICE_BCC:
            return 14;
        case LATTICE_CUBIC_DIAMOND:
        case LATTICE_HEX_DIAMOND:
            return 16;
        case LATTICE_SC:
            return 6;
        default:
            return 0;
    }
}

inline StructureType structureTypeFor(CoordinationStructureType coordinationType){
    switch(coordinationType){
        case COORD_CUBIC_DIAMOND:
            return StructureType::CUBIC_DIAMOND;
        case COORD_HEX_DIAMOND:
            return StructureType::HEX_DIAMOND;
        case COORD_FCC:
            return StructureType::FCC;
        case COORD_HCP:
            return StructureType::HCP;
        case COORD_BCC:
            return StructureType::BCC;
        case COORD_SC:
            return StructureType::SC;
        default:
            return StructureType::OTHER;
    }
}

inline bool buildExpandedNeighborShell(
    const NativeCnaLocalEnvironmentDescriptor& descriptor,
    const NearestNeighborFinder& neighList,
    const NearestNeighborFinder::Query<MAX_NEIGHBORS>& neighQuery,
    int particleIndex,
    int coordinationNumber,
    int* neighborIndices,
    Vector3* neighborVectors,
    NeighborBondArray& neighborArray
){
    if(descriptor.expansionSeedCount <= 0 || coordinationNumber <= descriptor.expansionSeedCount){
        return false;
    }

    const int neighborsPerSeed =
        (coordinationNumber - descriptor.expansionSeedCount) / descriptor.expansionSeedCount;
    if(descriptor.expansionSeedCount * neighborsPerSeed + descriptor.expansionSeedCount != coordinationNumber){
        return false;
    }

    int outputIndex = descriptor.expansionSeedCount;
    for(int seedIndex = 0; seedIndex < descriptor.expansionSeedCount; ++seedIndex){
        const Vector3& seedVector = neighQuery.results()[seedIndex].delta;
        neighborVectors[seedIndex] = seedVector;
        neighborIndices[seedIndex] = neighQuery.results()[seedIndex].index;

        NearestNeighborFinder::Query<MAX_NEIGHBORS> seedQuery(neighList);
        seedQuery.findNeighbors(neighList.particlePos(neighborIndices[seedIndex]));
        if(static_cast<int>(seedQuery.results().size()) < descriptor.expansionSeedCount){
            return false;
        }

        int produced = 0;
        for(int neighborIndex = 0; neighborIndex < descriptor.expansionSeedCount; ++neighborIndex){
            Vector3 vector = seedVector + seedQuery.results()[neighborIndex].delta;
            if(seedQuery.results()[neighborIndex].index == particleIndex && vector.isZero()){
                continue;
            }
            if(outputIndex == coordinationNumber){
                return false;
            }

            neighborIndices[outputIndex] = seedQuery.results()[neighborIndex].index;
            neighborVectors[outputIndex] = vector;
            neighborArray.setNeighborBond(seedIndex, outputIndex, true);
            ++outputIndex;
            ++produced;
        }

        if(produced != neighborsPerSeed){
            return false;
        }
    }

    return outputIndex == coordinationNumber;
}

inline double computeLocalCutoff(
    LatticeStructureType inputCrystalType,
    const NearestNeighborFinder& neighList,
    const NearestNeighborFinder::Query<MAX_NEIGHBORS>& neighQuery,
    int numNeighbors,
    int coordinationNumber,
    int particleIndex,
    int* neighborIndices,
    Vector3* neighborVectors,
    NeighborBondArray& neighborArray,
    bool allowCollapsedNeighborShell = false,
    bool* outRejectedByExtraNeighbor = nullptr
){
    if(outRejectedByExtraNeighbor){
        *outRejectedByExtraNeighbor = false;
    }

    const NativeCnaLocalEnvironmentDescriptor descriptor = nativeCnaLocalEnvironmentFor(inputCrystalType);
    if(descriptor.construction == NativeCnaLocalEnvironmentConstruction::None ||
       descriptor.referenceNeighborCount <= 0 ||
       descriptor.cutoffMultiplier <= 0.0){
        return 0.0;
    }

    if(descriptor.construction == NativeCnaLocalEnvironmentConstruction::Direct){
        if(numNeighbors < coordinationNumber || numNeighbors < descriptor.referenceNeighborCount){
            return 0.0;
        }

        for(int neighborIndex = 0; neighborIndex < coordinationNumber; ++neighborIndex){
            neighborIndices[neighborIndex] = neighQuery.results()[neighborIndex].index;
            neighborVectors[neighborIndex] = neighQuery.results()[neighborIndex].delta;
        }
    }else if(descriptor.construction == NativeCnaLocalEnvironmentConstruction::FirstShellExpansion){
        if(numNeighbors < descriptor.expansionSeedCount){
            return 0.0;
        }
        if(!buildExpandedNeighborShell(
            descriptor,
            neighList,
            neighQuery,
            particleIndex,
            coordinationNumber,
            neighborIndices,
            neighborVectors,
            neighborArray
        )){
            return 0.0;
        }
    }else{
        return 0.0;
    }

    const int referenceBegin = descriptor.referenceNeighborOffset;
    const int referenceEnd = referenceBegin + descriptor.referenceNeighborCount;
    if(referenceBegin < 0 || referenceEnd > coordinationNumber){
        return 0.0;
    }

    double localScaling = 0.0;
    for(int neighborIndex = referenceBegin; neighborIndex < referenceEnd; ++neighborIndex){
        localScaling += neighborVectors[neighborIndex].length();
    }
    localScaling /= static_cast<double>(descriptor.referenceNeighborCount);

    const double localCutoff = localScaling * descriptor.cutoffMultiplier;
    const double localCutoffSquared = localCutoff * localCutoff;

    if(descriptor.extraNeighborRejectIndex >= 0 &&
       numNeighbors > descriptor.extraNeighborRejectIndex &&
       neighQuery.results()[descriptor.extraNeighborRejectIndex].distanceSq <= localCutoffSquared){
        if(outRejectedByExtraNeighbor){
            *outRejectedByExtraNeighbor = true;
        }
        if(!allowCollapsedNeighborShell){
            return 0.0;
        }
    }

    if(descriptor.construction == NativeCnaLocalEnvironmentConstruction::Direct){
        for(int ni1 = 0; ni1 < coordinationNumber; ++ni1){
            neighborArray.setNeighborBond(ni1, ni1, false);
            for(int ni2 = ni1 + 1; ni2 < coordinationNumber; ++ni2){
                neighborArray.setNeighborBond(
                    ni1,
                    ni2,
                    (neighborVectors[ni1] - neighborVectors[ni2]).squaredLength() <= localCutoffSquared
                );
            }
        }
    }else{
        for(int ni1 = descriptor.bondStartIndex; ni1 < coordinationNumber; ++ni1){
            for(int ni2 = ni1 + 1; ni2 < coordinationNumber; ++ni2){
                const Vector3 distance = neighborVectors[ni1] - neighborVectors[ni2];
                neighborArray.setNeighborBond(ni1, ni2, distance.squaredLength() <= localCutoffSquared);
            }
        }
    }

    return localCutoff;
}

inline bool determineLocalStructure(
    LatticeStructureType inputCrystalType,
    bool identifyPlanarDefects,
    const NearestNeighborFinder& neighList,
    int particleIndex,
    const CoordinationStructure* coordinationStructures,
    LocalStructureMatch& outMatch
){
    outMatch = {};

    std::array<int, MAX_NEIGHBORS> cnaSignatures;
    std::array<int, MAX_NEIGHBORS> previousMapping;
    NeighborBondArray neighborArray;

    NearestNeighborFinder::Query<MAX_NEIGHBORS> neighQuery(neighList);
    neighQuery.findNeighbors(neighList.particlePos(particleIndex));
    const int numNeighbors = static_cast<int>(neighQuery.results().size());

    outMatch.coordinationNumber = coordinationNumberFor(inputCrystalType);
    if(outMatch.coordinationNumber <= 0 || numNeighbors < outMatch.coordinationNumber){
        return false;
    }

    bool rejectedByExtraNeighbor = false;
    outMatch.localCutoff = computeLocalCutoff(
        inputCrystalType,
        neighList,
        neighQuery,
        numNeighbors,
        outMatch.coordinationNumber,
        particleIndex,
        outMatch.neighborIndices.data(),
        outMatch.neighborVectors.data(),
        neighborArray,
        false,
        &rejectedByExtraNeighbor
    );
    if(outMatch.localCutoff == 0.0 && rejectedByExtraNeighbor){
        outMatch.localCutoff = computeLocalCutoff(
            inputCrystalType,
            neighList,
            neighQuery,
            numNeighbors,
            outMatch.coordinationNumber,
            particleIndex,
            outMatch.neighborIndices.data(),
            outMatch.neighborVectors.data(),
            neighborArray,
            true,
            nullptr
        );
    }
    if(outMatch.localCutoff == 0.0){
        return false;
    }

    for(int n = 0; n < outMatch.coordinationNumber; ++n){
        outMatch.neighborMapping[n] = n;
        previousMapping[n] = -1;
    }

    outMatch.coordinationType = CommonNeighborAnalysis::computeCoordinationType(
        neighborArray,
        outMatch.coordinationNumber,
        cnaSignatures.data(),
        inputCrystalType,
        identifyPlanarDefects
    );
    if(outMatch.coordinationType == COORD_OTHER){
        return false;
    }

    const bool found = CommonNeighborAnalysis::findMatchingNeighborPermutation(
        outMatch.coordinationType,
        outMatch.neighborMapping.data(),
        previousMapping.data(),
        outMatch.coordinationNumber,
        cnaSignatures.data(),
        neighborArray,
        coordinationStructures
    );
    if(!found){
        return false;
    }

    outMatch.structureType = structureTypeFor(outMatch.coordinationType);
    return outMatch.structureType != StructureType::OTHER;
}

}
