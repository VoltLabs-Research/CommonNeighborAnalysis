#include <volt/analysis/structure_analysis.h>
#include <volt/coordination_structures.h>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

#include <algorithm>
#include <limits>
#include <mutex>
#include <vector>

namespace Volt {

namespace {

void ensureCoordinationStructuresInitialized(){
    static std::once_flag initFlag;
    std::call_once(initFlag, []() {
        CoordinationStructures::initializeStructures();
    });
}

}

int StructureAnalysis::findClosestSymmetryPermutation(int structureType, const Matrix3& rotation){
    ensureCoordinationStructuresInitialized();

    const LatticeStructure& lattice = CoordinationStructures::getLatticeStruct(structureType);
    int bestIndex = 0;
    double bestDeviation = std::numeric_limits<double>::max();

    for(int i = 0; i < lattice.permutations.size(); ++i){
        const Matrix3& sym = lattice.permutations[i].transformation;
        double deviation = 0;
        for(int r = 0; r < 3; ++r){
            for(int c = 0; c < 3; ++c){
                double diff = rotation(r, c) - sym(r, c);
                deviation += diff * diff;
            }
        }
        if(deviation < bestDeviation){
            bestDeviation = deviation;
            bestIndex = i;
        }
    }
    return bestIndex;
}

int StructureAnalysis::coordinationNumber(int structureType) const{
    ensureCoordinationStructuresInitialized();
    return CoordinationStructures::getCoordStruct(structureType).numNeighbors;
}

int StructureAnalysis::commonNeighborIndex(int structureType, int neighborIndex, int commonNeighborSlot) const{
    ensureCoordinationStructuresInitialized();
    return CoordinationStructures::getCoordStruct(structureType).commonNeighbors[neighborIndex][commonNeighborSlot];
}

int StructureAnalysis::symmetryPermutationCount(int structureType) const{
    ensureCoordinationStructuresInitialized();
    return static_cast<int>(CoordinationStructures::getLatticeStruct(structureType).permutations.size());
}

int StructureAnalysis::symmetryPermutationEntry(int structureType, int symmetryIndex, int neighborIndex) const{
    ensureCoordinationStructuresInitialized();
    return CoordinationStructures::getLatticeStruct(structureType).permutations[symmetryIndex].permutation[neighborIndex];
}

const Matrix3& StructureAnalysis::symmetryTransformation(int structureType, int symmetryIndex) const{
    ensureCoordinationStructuresInitialized();
    return CoordinationStructures::getLatticeStruct(structureType).permutations[symmetryIndex].transformation;
}

int StructureAnalysis::symmetryInverseProduct(int structureType, int symmetryIndex, int transformationIndex) const{
    ensureCoordinationStructuresInitialized();
    return CoordinationStructures::getLatticeStruct(structureType).permutations[symmetryIndex].inverseProduct[transformationIndex];
}

const Vector3& StructureAnalysis::latticeVector(int structureType, int latticeVectorIndex) const{
    ensureCoordinationStructuresInitialized();
    return CoordinationStructures::getLatticeStruct(structureType).latticeVectors[latticeVectorIndex];
}

const Vector3& StructureAnalysis::neighborLatticeVector(int centralAtomIndex, int neighborIndex) const{
    ensureCoordinationStructuresInitialized();

    assert(_context.atomSymmetryPermutations);
    const int structureType = _context.structureTypes->getInt(centralAtomIndex);
    assert(neighborIndex >= 0 && neighborIndex < coordinationNumber(structureType));
    const int symmetryPermutationIndex = _context.atomSymmetryPermutations->getInt(centralAtomIndex);
    assert(symmetryPermutationIndex >= 0 && symmetryPermutationIndex < symmetryPermutationCount(structureType));
    return latticeVector(structureType, symmetryPermutationEntry(structureType, symmetryPermutationIndex, neighborIndex));
}

void StructureAnalysis::identifyStructuresCNA(){
    ensureCoordinationStructuresInitialized();

    const int maxNeighborListSize = MAX_NEIGHBORS;
    NearestNeighborFinder neighFinder(maxNeighborListSize);
    if(!neighFinder.prepare(_context.positions, _context.simCell, _context.particleSelection)){
        throw std::runtime_error("Error in neighFinder.preapre(...)");
    }

    CoordinationStructures coordinationStructures(
        _context.structureTypes,
        _context.inputCrystalType,
        _identifyPlanarDefects,
        _context.simCell
    );

    const size_t N = _context.atomCount();
    std::vector<int> localCounts(N, 0);
    std::vector<std::array<int, MAX_NEIGHBORS>> orderedNeighborIndices(N);
    _maximumNeighborDistance = tbb::parallel_reduce(
        tbb::blocked_range<size_t>(0, N),
        0.0,
        [this, &neighFinder, &coordinationStructures, &localCounts, &orderedNeighborIndices](const tbb::blocked_range<size_t>& range, double maxDistSoFar) -> double {
            for(size_t index = range.begin(); index != range.end(); ++index){
                int count = 0;
                auto& neighbors = orderedNeighborIndices[index];
                neighbors.fill(-1);
                double localMaxDistance = coordinationStructures.determineLocalStructure(
                    neighFinder,
                    index,
                    &count,
                    neighbors.data()
                );
                localCounts[index] = count;
                if(localMaxDistance > maxDistSoFar){
                    maxDistSoFar = localMaxDistance;
                }
            }
            return maxDistSoFar;
        },
        [](double a, double b) -> double {
            return std::max(a, b);
        }
    );

    auto* offsets = _context.neighborOffsets->dataInt();
    offsets[0] = 0;
    for(size_t i = 0; i < N; ++i){
        offsets[i + 1] = offsets[i] + localCounts[i];
    }
    const size_t totalNeighbors = static_cast<size_t>(offsets[N]);
    _context.neighborIndices = std::make_shared<ParticleProperty>(
        totalNeighbors, DataType::Int, 1, 0, false);

    auto* indices = _context.neighborIndices->dataInt();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, N), [&](const auto& range){
        for(size_t index = range.begin(); index != range.end(); ++index){
            const int count = localCounts[index];
            if(count <= 0){
                continue;
            }
            const int start = offsets[index];
            for(int j = 0; j < count; ++j){
                indices[start + j] = orderedNeighborIndices[index][j];
            }
            _context.neighborCounts->setInt(index, count);
        }
    });
}

}
