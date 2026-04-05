#pragma once

#include <volt/analysis/nearest_neighbor_finder.h>
#include <volt/analysis/cna_classifier.h>
#include <volt/topology/crystal_coordination_pattern.h>

#include <array>

namespace Volt::CnaLocalStructureUtils{

struct LocalStructureMatch{
    double localCutoff = 0.0;
    int coordinationNumber = 0;
    CoordinationStructureType coordinationType = COORD_OTHER;
    StructureType structureType = StructureType::OTHER;
    std::array<int, MAX_NEIGHBORS> neighborIndices{};
    std::array<int, MAX_NEIGHBORS> neighborMapping{};
};

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

inline double computeLocalCutoff(
    LatticeStructureType inputCrystalType,
    const NearestNeighborFinder& neighList,
    const NearestNeighborFinder::Query<MAX_NEIGHBORS>& neighQuery,
    int numNeighbors,
    int coordinationNumber,
    int particleIndex,
    int* neighborIndices,
    Vector3* neighborVectors,
    NeighborBondArray& neighborArray
){
    double localScaling = 0.0;
    double localCutoff = 0.0;

    switch(inputCrystalType){
        case LATTICE_FCC:
        case LATTICE_HCP:
            for(int neighbor = 0; neighbor < 12; ++neighbor){
                localScaling += sqrt(neighQuery.results()[neighbor].distanceSq);
            }
            localScaling /= 12;
            localCutoff = localScaling * (1.0f + sqrt(2.0f)) * 0.5f;
            break;
        case LATTICE_BCC:
            for(int neighbor = 0; neighbor < 8; ++neighbor){
                localScaling += sqrt(neighQuery.results()[neighbor].distanceSq);
            }
            localScaling /= 8;
            localCutoff = localScaling / (sqrt(3.0) / 2.0) * 0.5 * (1.0 + sqrt(2.0));
            break;
        case LATTICE_CUBIC_DIAMOND:
        case LATTICE_HEX_DIAMOND: {
            int outputIndex = 4;
            for(int i = 0; i < 4; ++i){
                const Vector3& v0 = neighQuery.results()[i].delta;
                neighborVectors[i] = v0;
                neighborIndices[i] = neighQuery.results()[i].index;

                NearestNeighborFinder::Query<MAX_NEIGHBORS> neighQuery2(neighList);
                neighQuery2.findNeighbors(neighList.particlePos(neighborIndices[i]));
                if(neighQuery2.results().size() < 4){
                    return 0.0;
                }

                for(int j = 0; j < 4; ++j){
                    Vector3 v = v0 + neighQuery2.results()[j].delta;
                    if(neighQuery2.results()[j].index == particleIndex && v.isZero()){
                        continue;
                    }
                    if(outputIndex == 16){
                        return 0.0;
                    }

                    neighborIndices[outputIndex] = neighQuery2.results()[j].index;
                    neighborVectors[outputIndex] = v;
                    neighborArray.setNeighborBond(i, outputIndex, true);
                    ++outputIndex;
                }

                if(outputIndex != (i * 3) + 7){
                    return 0.0;
                }
            }

            for(int neighbor = 4; neighbor < 16; ++neighbor){
                localScaling += neighborVectors[neighbor].length();
            }

            localScaling /= 12;
            localCutoff = localScaling * 1.2071068;
            break;
        }
        default:
            return 0.0;
    }

    const double localCutoffSquared = localCutoff * localCutoff;

    switch(inputCrystalType){
        case LATTICE_FCC:
        case LATTICE_HCP:
        case LATTICE_BCC:
        case LATTICE_SC:
            if(numNeighbors > coordinationNumber && neighQuery.results()[coordinationNumber].distanceSq <= localCutoffSquared){
                return 0.0;
            }

            for(int ni1 = 0; ni1 < coordinationNumber; ++ni1){
                neighborIndices[ni1] = neighQuery.results()[ni1].index;
                neighborVectors[ni1] = neighQuery.results()[ni1].delta;
                neighborArray.setNeighborBond(ni1, ni1, false);
                for(int ni2 = ni1 + 1; ni2 < coordinationNumber; ++ni2){
                    neighborArray.setNeighborBond(
                        ni1,
                        ni2,
                        (neighQuery.results()[ni1].delta - neighQuery.results()[ni2].delta).squaredLength() <= localCutoffSquared
                    );
                }
            }
            break;
        case LATTICE_CUBIC_DIAMOND:
        case LATTICE_HEX_DIAMOND:
            for(int ni1 = 4; ni1 < coordinationNumber; ++ni1){
                for(int ni2 = ni1 + 1; ni2 < coordinationNumber; ++ni2){
                    const Vector3 distance = neighborVectors[ni1] - neighborVectors[ni2];
                    neighborArray.setNeighborBond(ni1, ni2, distance.squaredLength() <= localCutoffSquared);
                }
            }
            break;
        default:
            return 0.0;
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

    std::array<Vector3, MAX_NEIGHBORS> neighborVectors;
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

    outMatch.localCutoff = computeLocalCutoff(
        inputCrystalType,
        neighList,
        neighQuery,
        numNeighbors,
        outMatch.coordinationNumber,
        particleIndex,
        outMatch.neighborIndices.data(),
        neighborVectors.data(),
        neighborArray
    );
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
