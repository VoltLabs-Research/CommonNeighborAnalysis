#include <volt/analysis/cna_structure_analysis.h>
#include <volt/analysis/crystal_topology_library.h>
#include <volt/topology/crystal_coordination_topology_init.h>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

#include <algorithm>
#include <vector>

namespace Volt {

namespace CnaStructureAnalysisDetail{

class CnaCrystalInfoProvider final : public StructureAnalysisCrystalInfo{
public:
    int findClosestSymmetryPermutation(int structureType, const Matrix3& rotation) const override{
        const SharedCrystalTopology* topology = sharedCrystalTopology(structureType);
        return topology ? findClosestSharedCrystalSymmetryPermutation(*topology, rotation) : 0;
    }

    int coordinationNumber(int structureType) const override{
        const SharedCrystalTopology* topology = sharedCrystalTopology(structureType);
        return topology ? topology->coordinationNumber : 0;
    }

    int commonNeighborIndex(int structureType, int neighborIndex, int commonNeighborSlot) const override{
        const SharedCrystalTopology* topology = sharedCrystalTopology(structureType);
        if(!topology ||
           neighborIndex < 0 ||
           neighborIndex >= topology->coordinationNumber ||
           commonNeighborSlot < 0 ||
           commonNeighborSlot > 1){
            return -1;
        }
        return topology->commonNeighbors[static_cast<std::size_t>(neighborIndex)]
            [static_cast<std::size_t>(commonNeighborSlot)];
    }

    int symmetryPermutationCount(int structureType) const override{
        const SharedCrystalTopology* topology = sharedCrystalTopology(structureType);
        return topology ? static_cast<int>(topology->symmetries.size()) : 0;
    }

    int symmetryPermutationEntry(int structureType, int symmetryIndex, int neighborIndex) const override{
        const SharedCrystalTopology* topology = sharedCrystalTopology(structureType);
        if(!topology ||
           symmetryIndex < 0 ||
           symmetryIndex >= static_cast<int>(topology->symmetries.size()) ||
           neighborIndex < 0 ||
           neighborIndex >= topology->coordinationNumber){
            return neighborIndex;
        }
        return topology->symmetries[static_cast<std::size_t>(symmetryIndex)]
            .permutation[static_cast<std::size_t>(neighborIndex)];
    }

    const Matrix3& symmetryTransformation(int structureType, int symmetryIndex) const override{
        const SharedCrystalTopology* topology = sharedCrystalTopology(structureType);
        static const Matrix3 identity = Matrix3::Identity();
        if(!topology ||
           symmetryIndex < 0 ||
           symmetryIndex >= static_cast<int>(topology->symmetries.size())){
            return identity;
        }
        return topology->symmetries[static_cast<std::size_t>(symmetryIndex)].transformation;
    }

    int symmetryInverseProduct(int structureType, int symmetryIndex, int transformationIndex) const override{
        const SharedCrystalTopology* topology = sharedCrystalTopology(structureType);
        if(!topology ||
           symmetryIndex < 0 ||
           symmetryIndex >= static_cast<int>(topology->symmetries.size())){
            return 0;
        }
        const auto& inverseProduct = topology->symmetries[static_cast<std::size_t>(symmetryIndex)].inverseProduct;
        if(transformationIndex < 0 || transformationIndex >= static_cast<int>(inverseProduct.size())){
            return 0;
        }
        return inverseProduct[static_cast<std::size_t>(transformationIndex)];
    }

    const Vector3& latticeVector(int structureType, int latticeVectorIndex) const override{
        static const Vector3 zero = Vector3::Zero();
        const SharedCrystalTopology* topology = sharedCrystalTopology(structureType);
        if(!topology ||
           latticeVectorIndex < 0 ||
           latticeVectorIndex >= static_cast<int>(topology->latticeVectors.size())){
            return zero;
        }
        return topology->latticeVectors[static_cast<std::size_t>(latticeVectorIndex)];
    }
};

std::shared_ptr<const StructureAnalysisCrystalInfo> cnaCrystalInfoProvider(){
    static const auto provider = std::make_shared<CnaCrystalInfoProvider>();
    return provider;
}

}

using namespace CnaStructureAnalysisDetail;

void identifyStructuresCNA(StructureAnalysis& analysis){
    StructureContext& context = analysis.context();
    ensureCoordinationStructuresInitialized();
    analysis.setCrystalInfoProvider(cnaCrystalInfoProvider());

    const int maxNeighborListSize = MAX_NEIGHBORS;
    NearestNeighborFinder neighFinder(maxNeighborListSize);
    if(!neighFinder.prepare(context.positions, context.simCell, context.particleSelection)){
        throw std::runtime_error("Error in neighFinder.preapre(...)");
    }

    CoordinationStructures coordinationStructures(
        context.structureTypes,
        context.inputCrystalType,
        true,
        context.simCell
    );

    const size_t N = context.atomCount();
    std::vector<int> localCounts(N, 0);
    std::vector<std::array<int, MAX_NEIGHBORS>> orderedNeighborIndices(N);
    context.maximumNeighborDistance = tbb::parallel_reduce(
        tbb::blocked_range<size_t>(0, N),
        0.0,
        [&neighFinder, &coordinationStructures, &localCounts, &orderedNeighborIndices](const tbb::blocked_range<size_t>& range, double maxDistSoFar) -> double {
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

    auto* offsets = context.neighborOffsets->dataInt();
    offsets[0] = 0;
    for(size_t i = 0; i < N; ++i){
        offsets[i + 1] = offsets[i] + localCounts[i];
    }
    const size_t totalNeighbors = static_cast<size_t>(offsets[N]);
    context.neighborIndices = std::make_shared<ParticleProperty>(
        totalNeighbors, DataType::Int, 1, 0, false);

    auto* indices = context.neighborIndices->dataInt();
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
            context.neighborCounts->setInt(index, count);
        }
    });

    std::vector<Vector3> neighborVectorOverrides(
        N * static_cast<std::size_t>(MAX_NEIGHBORS),
        Vector3::Zero()
    );
    for(std::size_t atomIndex = 0; atomIndex < N; ++atomIndex){
        const int structureType = context.structureTypes->getInt(atomIndex);
        if(structureType == LATTICE_OTHER){
            continue;
        }

        const SharedCrystalTopology* topology = sharedCrystalTopology(structureType);
        const CrystalTopologyEntry* topologyEntry = crystalTopologyByStructureType(structureType);
        if(!topology || !topologyEntry){
            continue;
        }

        const int exportableCount = std::min(localCounts[atomIndex], topology->coordinationNumber);
        const int exportSymmetryIndex = topologyEntry->exportSymmetryIndex >= 0 &&
            topologyEntry->exportSymmetryIndex < static_cast<int>(topology->symmetries.size())
            ? topologyEntry->exportSymmetryIndex
            : 0;
        for(int neighborSlot = 0; neighborSlot < exportableCount; ++neighborSlot){
            int latticeVectorIndex = neighborSlot;
            if(!topology->symmetries.empty()){
                latticeVectorIndex = topology->symmetries[static_cast<std::size_t>(exportSymmetryIndex)]
                    .permutation[static_cast<std::size_t>(neighborSlot)];
            }
            if(latticeVectorIndex < 0 || latticeVectorIndex >= static_cast<int>(topology->latticeVectors.size())){
                continue;
            }
            neighborVectorOverrides[
                atomIndex * static_cast<std::size_t>(MAX_NEIGHBORS) +
                static_cast<std::size_t>(neighborSlot)
            ] = topology->latticeVectors[static_cast<std::size_t>(latticeVectorIndex)];
        }
    }
    analysis.setNeighborLatticeVectorOverrides(
        std::move(neighborVectorOverrides),
        static_cast<std::size_t>(MAX_NEIGHBORS)
    );
}

}
