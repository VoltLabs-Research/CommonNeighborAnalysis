#include <volt/cna_cluster_input_adapter.h>

#include <algorithm>
#include <cstdint>
#include <limits>

namespace Volt{

namespace{

std::uint64_t fullSymmetryMask(int symmetryCount){
    if(symmetryCount <= 0){
        return 0;
    }
    if(symmetryCount >= 63){
        return std::numeric_limits<std::uint64_t>::max();
    }
    return (std::uint64_t{1} << symmetryCount) - 1;
}

}

void CNAClusterInputAdapter::prepare(StructureAnalysis& analysis, AnalysisContext& context){
    analysis.setClusterRuleProvider(nullptr);

    if(!context.atomAllowedSymmetryMasks){
        context.atomAllowedSymmetryMasks = std::make_shared<ParticleProperty>(
            context.atomCount(),
            DataType::Int64,
            1,
            0,
            true
        );
    }

    std::fill(
        context.atomSymmetryPermutations->dataInt(),
        context.atomSymmetryPermutations->dataInt() + context.atomSymmetryPermutations->size(),
        -1
    );
    std::fill(
        context.atomAllowedSymmetryMasks->dataInt64(),
        context.atomAllowedSymmetryMasks->dataInt64() + context.atomAllowedSymmetryMasks->size(),
        0
    );

    for(size_t atomIndex = 0; atomIndex < context.atomCount(); ++atomIndex){
        const int structureType = context.structureTypes->getInt(atomIndex);
        if(structureType == LATTICE_OTHER){
            continue;
        }

        context.atomAllowedSymmetryMasks->setInt64(
            atomIndex,
            static_cast<std::int64_t>(fullSymmetryMask(analysis.symmetryPermutationCount(structureType)))
        );
    }
}

}
