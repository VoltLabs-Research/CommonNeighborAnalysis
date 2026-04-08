#pragma once

#include <volt/analysis/crystal_topology_library.h>
#include <volt/analysis/cna_classifier.h>
#include <volt/analysis/nearest_neighbor_finder.h>
#include <volt/math/quaternion.h>
#include <volt/structures/crystal_topology_registry.h>
#include <volt/topology/crystal_coordination_pattern.h>

#include <algorithm>
#include <array>
#include <cstdlib>
#include <limits>
#include <spdlog/spdlog.h>

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
    const auto* topology = crystalTopologyByLatticeType(static_cast<int>(inputCrystalType));
    return topology ? topology->coordinationNumber : 0;
}

inline StructureType structureTypeFor(CoordinationStructureType coordinationType){
    const auto* topology = crystalTopologyByCoordinationType(static_cast<int>(coordinationType));
    return topology ? static_cast<StructureType>(topology->structureType) : StructureType::OTHER;
}

inline bool cnaLocalDebugEnabled(){
    static const bool enabled = []() {
        const char* value = std::getenv("VOLT_DEBUG_CNA_LOCAL");
        return value != nullptr && value[0] != '\0' && value[0] != '0';
    }();
    return enabled;
}

inline bool orthonormalizeOrientation(const Matrix3& input, Matrix3& output){
    Vector3 c0 = input.column(0);
    Vector3 c1 = input.column(1);
    Vector3 c2 = input.column(2);

    const double l0 = c0.length();
    if(l0 <= EPSILON){
        return false;
    }
    c0 /= l0;

    c1 -= c0 * c0.dot(c1);
    const double l1 = c1.length();
    if(l1 <= EPSILON){
        return false;
    }
    c1 /= l1;

    c2 -= c0 * c0.dot(c2);
    c2 -= c1 * c1.dot(c2);
    const double l2 = c2.length();
    if(l2 <= EPSILON){
        return false;
    }
    c2 /= l2;

    output = Matrix3(c0, c1, c2);
    if(output.determinant() < 0.0){
        output.column(2) = -output.column(2);
    }
    return true;
}

inline bool topologySupportsOrientationInversion(
    const SharedCrystalTopology& topology,
    std::array<int, MAX_NEIGHBORS>& outInversePermutation
){
    outInversePermutation.fill(-1);
    for(int slot = 0; slot < topology.coordinationNumber; ++slot){
        const Vector3& vector = topology.latticeVectors[static_cast<std::size_t>(slot)];
        for(int candidateSlot = 0; candidateSlot < topology.coordinationNumber; ++candidateSlot){
            if((topology.latticeVectors[static_cast<std::size_t>(candidateSlot)] + vector).isZero(EPSILON)){
                outInversePermutation[static_cast<std::size_t>(slot)] = candidateSlot;
                break;
            }
        }
        if(outInversePermutation[static_cast<std::size_t>(slot)] < 0){
            return false;
        }
    }
    return true;
}

inline bool computeMatchedOrientation(
    const SharedCrystalTopology& topology,
    const Vector3* neighborVectors,
    const int* neighborMapping,
    Matrix3& outRawOrientation,
    Matrix3& outOrientation
){
    Matrix3 orientationV = Matrix3::Zero();
    Matrix3 orientationW = Matrix3::Zero();
    int vectorCount = 0;

    for(int canonicalSlot = 0; canonicalSlot < topology.coordinationNumber; ++canonicalSlot){
        const int localSlot = neighborMapping[canonicalSlot];
        if(localSlot < 0 || localSlot >= topology.coordinationNumber){
            return false;
        }

        const Vector3& idealVector = topology.latticeVectors[static_cast<std::size_t>(canonicalSlot)];
        const Vector3& spatialVector = neighborVectors[localSlot];
        for(int row = 0; row < 3; ++row){
            for(int column = 0; column < 3; ++column){
                orientationV(row, column) += idealVector[column] * idealVector[row];
                orientationW(row, column) += idealVector[column] * spatialVector[row];
            }
        }
        ++vectorCount;
    }

    if(vectorCount < 3){
        return false;
    }

    Matrix3 orientationVInverse;
    if(!orientationV.inverse(orientationVInverse)){
        return false;
    }

    outRawOrientation = Matrix3(orientationW * orientationVInverse);
    return orthonormalizeOrientation(outRawOrientation, outOrientation);
}

inline bool computeMatchedOrientation(
    const SharedCrystalTopology& topology,
    const Vector3* neighborVectors,
    const int* neighborMapping,
    Matrix3& outOrientation
){
    Matrix3 rawOrientation;
    return computeMatchedOrientation(
        topology,
        neighborVectors,
        neighborMapping,
        rawOrientation,
        outOrientation
    );
}

inline double identityDeviation(const Matrix3& matrix){
    double error = 0.0;
    for(int row = 0; row < 3; ++row){
        for(int column = 0; column < 3; ++column){
            const double expected = row == column ? 1.0 : 0.0;
            const double delta = matrix(row, column) - expected;
            error += delta * delta;
        }
    }
    return error;
}

inline bool canonicalizeNeighborMapping(
    StructureType structureType,
    const Vector3* neighborVectors,
    int coordinationNumber,
    int* neighborMapping
){
    const SharedCrystalTopology* topology = sharedCrystalTopology(static_cast<int>(structureType));
    if(!topology || coordinationNumber != topology->coordinationNumber){
        return true;
    }

    auto scoreOrientation = [](const Matrix3& orientation) {
        Quaternion q(orientation);
        q.normalize();
        if(q.w() < 0.0){
            q = -q;
        }
        return std::array<double, 4>{q.w(), q.x(), q.y(), q.z()};
    };

    auto isBetterScore = [](const std::array<double, 4>& lhs, const std::array<double, 4>& rhs) {
        constexpr double tolerance = 1e-12;
        for(std::size_t index = 0; index < lhs.size(); ++index){
            if(lhs[index] > rhs[index] + tolerance){
                return true;
            }
            if(lhs[index] + tolerance < rhs[index]){
                return false;
            }
        }
        return false;
    };

    std::array<int, MAX_NEIGHBORS> baseMapping{};
    std::copy_n(neighborMapping, coordinationNumber, baseMapping.begin());

    std::array<int, MAX_NEIGHBORS> inversePermutation{};
    const bool allowInversionEquivalent = topologySupportsOrientationInversion(
        *topology,
        inversePermutation
    );

    std::array<int, MAX_NEIGHBORS> bestMapping{};
    std::array<double, 4> bestScore{};
    bool bestValid = false;
    bool bestProper = false;

    auto tryCandidate = [&](const std::array<int, MAX_NEIGHBORS>& canonicalPermutation) {
        std::array<int, MAX_NEIGHBORS> candidateMapping{};
        for(int canonicalSlot = 0; canonicalSlot < coordinationNumber; ++canonicalSlot){
            const int mappedCanonicalSlot = canonicalPermutation[static_cast<std::size_t>(canonicalSlot)];
            if(mappedCanonicalSlot < 0 || mappedCanonicalSlot >= coordinationNumber){
                return;
            }
            candidateMapping[static_cast<std::size_t>(canonicalSlot)] =
                baseMapping[static_cast<std::size_t>(mappedCanonicalSlot)];
        }

        Matrix3 rawOrientation;
        Matrix3 orientation;
        if(!computeMatchedOrientation(
            *topology,
            neighborVectors,
            candidateMapping.data(),
            rawOrientation,
            orientation
        )){
            return;
        }

        const bool properOrientation = rawOrientation.determinant() >= 0.0;
        const std::array<double, 4> score = scoreOrientation(orientation);
        if(!bestValid ||
           (properOrientation && !bestProper) ||
           (properOrientation == bestProper && isBetterScore(score, bestScore))){
            bestMapping = candidateMapping;
            bestScore = score;
            bestValid = true;
            bestProper = properOrientation;
        }
    };

    for(const auto& symmetry : topology->symmetries){
        tryCandidate(symmetry.permutation);
    }

    if(allowInversionEquivalent){
        for(const auto& symmetry : topology->symmetries){
            std::array<int, MAX_NEIGHBORS> invertedPermutation{};
            invertedPermutation.fill(-1);
            for(int canonicalSlot = 0; canonicalSlot < coordinationNumber; ++canonicalSlot){
                const int rotatedSlot = symmetry.permutation[static_cast<std::size_t>(canonicalSlot)];
                if(rotatedSlot < 0 || rotatedSlot >= coordinationNumber){
                    continue;
                }
                invertedPermutation[static_cast<std::size_t>(canonicalSlot)] =
                    inversePermutation[static_cast<std::size_t>(rotatedSlot)];
            }
            tryCandidate(invertedPermutation);
        }
    }

    if(!bestValid){
        return false;
    }

    std::copy_n(bestMapping.begin(), coordinationNumber, neighborMapping);
    return true;
}

inline bool buildExpandedNeighborShell(
    const CrystalTopologyEntry& topology,
    const NearestNeighborFinder& neighList,
    const NearestNeighborFinder::Query<MAX_NEIGHBORS>& neighQuery,
    int particleIndex,
    int* neighborIndices,
    Vector3* neighborVectors,
    NeighborBondArray& neighborArray
){
    const int seedCount = topology.cnaLocalEnvironment.expansionSeedCount;
    if(seedCount <= 0 || topology.coordinationNumber <= seedCount){
        return false;
    }

    const int neighborsPerSeed = (topology.coordinationNumber - seedCount) / seedCount;
    if(seedCount * neighborsPerSeed + seedCount != topology.coordinationNumber){
        return false;
    }

    int outputIndex = seedCount;
    for(int seedIndex = 0; seedIndex < seedCount; ++seedIndex){
        const Vector3& seedVector = neighQuery.results()[seedIndex].delta;
        neighborVectors[seedIndex] = seedVector;
        neighborIndices[seedIndex] = neighQuery.results()[seedIndex].index;

        NearestNeighborFinder::Query<MAX_NEIGHBORS> seedQuery(neighList);
        seedQuery.findNeighbors(neighList.particlePos(neighborIndices[seedIndex]));
        if(static_cast<int>(seedQuery.results().size()) < seedCount){
            return false;
        }

        int produced = 0;
        for(int neighborIndex = 0; neighborIndex < seedCount; ++neighborIndex){
            Vector3 vector = seedVector + seedQuery.results()[neighborIndex].delta;
            if(seedQuery.results()[neighborIndex].index == particleIndex && vector.isZero()){
                continue;
            }
            if(outputIndex == topology.coordinationNumber){
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

    return outputIndex == topology.coordinationNumber;
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

    const auto* topology = crystalTopologyByLatticeType(static_cast<int>(inputCrystalType));
    if(!topology){
        return 0.0;
    }

    const auto& localEnvironment = topology->cnaLocalEnvironment;
    if(localEnvironment.construction == CnaLocalEnvironmentConstruction::None ||
       localEnvironment.referenceNeighborCount <= 0 ||
       localEnvironment.cutoffMultiplier <= 0.0){
        return 0.0;
    }

    if(localEnvironment.construction == CnaLocalEnvironmentConstruction::Direct){
        if(numNeighbors < coordinationNumber || numNeighbors < localEnvironment.referenceNeighborCount){
            return 0.0;
        }

        for(int neighborIndex = 0; neighborIndex < coordinationNumber; ++neighborIndex){
            neighborIndices[neighborIndex] = neighQuery.results()[neighborIndex].index;
            neighborVectors[neighborIndex] = neighQuery.results()[neighborIndex].delta;
        }
    }else if(localEnvironment.construction == CnaLocalEnvironmentConstruction::FirstShellExpansion){
        if(numNeighbors < localEnvironment.expansionSeedCount){
            return 0.0;
        }

        if(!buildExpandedNeighborShell(
            *topology,
            neighList,
            neighQuery,
            particleIndex,
            neighborIndices,
            neighborVectors,
            neighborArray
        )){
            return 0.0;
        }
    }else{
        return 0.0;
    }

    const int referenceBegin = localEnvironment.referenceNeighborOffset;
    const int referenceEnd = referenceBegin + localEnvironment.referenceNeighborCount;
    if(referenceBegin < 0 || referenceEnd > coordinationNumber){
        return 0.0;
    }

    double localScaling = 0.0;
    for(int neighborIndex = referenceBegin; neighborIndex < referenceEnd; ++neighborIndex){
        localScaling += neighborVectors[neighborIndex].length();
    }
    localScaling /= static_cast<double>(localEnvironment.referenceNeighborCount);

    const double localCutoff = localScaling * localEnvironment.cutoffMultiplier;
    const double localCutoffSquared = localCutoff * localCutoff;

    if(localEnvironment.extraNeighborRejectIndex >= 0 &&
       numNeighbors > localEnvironment.extraNeighborRejectIndex &&
       neighQuery.results()[localEnvironment.extraNeighborRejectIndex].distanceSq <= localCutoffSquared){
        if(outRejectedByExtraNeighbor){
            *outRejectedByExtraNeighbor = true;
        }
        if(!allowCollapsedNeighborShell){
            return 0.0;
        }
    }

    if(localEnvironment.construction == CnaLocalEnvironmentConstruction::Direct){
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
        for(int ni1 = localEnvironment.bondStartIndex; ni1 < coordinationNumber; ++ni1){
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

    bool rejectedByExtraNeighbor = false;
    outMatch.localCutoff = computeLocalCutoff(
        inputCrystalType,
        neighList,
        neighQuery,
        numNeighbors,
        outMatch.coordinationNumber,
        particleIndex,
        outMatch.neighborIndices.data(),
        neighborVectors.data(),
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
            neighborVectors.data(),
            neighborArray,
            true,
            nullptr
        );
        if(cnaLocalDebugEnabled() && outMatch.localCutoff > 0.0 && particleIndex < 3){
            spdlog::info(
                "CNA local debug: particle={} lattice_type={} accepted with collapsed-shell fallback coord_num={} local_cutoff={}",
                particleIndex,
                inputCrystalType,
                outMatch.coordinationNumber,
                outMatch.localCutoff
            );
        }
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
        if(cnaLocalDebugEnabled() && particleIndex < 3){
            spdlog::info(
                "CNA local debug: particle={} lattice_type={} failed coordination match coord_num={} local_cutoff={}",
                particleIndex,
                inputCrystalType,
                outMatch.coordinationNumber,
                outMatch.localCutoff
            );
        }
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
        if(cnaLocalDebugEnabled() && particleIndex < 3){
            spdlog::info(
                "CNA local debug: particle={} lattice_type={} failed neighbor permutation coordination_type={}",
                particleIndex,
                inputCrystalType,
                outMatch.coordinationType
            );
        }
        return false;
    }

    outMatch.structureType = structureTypeFor(outMatch.coordinationType);
    if(outMatch.structureType == StructureType::OTHER){
        return false;
    }

    if(!canonicalizeNeighborMapping(
        outMatch.structureType,
        neighborVectors.data(),
        outMatch.coordinationNumber,
        outMatch.neighborMapping.data()
    )){
        return false;
    }

    return outMatch.structureType != StructureType::OTHER;
}

}
