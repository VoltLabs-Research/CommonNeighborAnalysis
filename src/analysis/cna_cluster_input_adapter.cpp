#include <volt/analysis/cna_cluster_input_adapter.h>
#include <volt/analysis/cluster_input_preparation.h>

namespace Volt{

void CNAClusterInputAdapter::prepare(StructureAnalysis& analysis, AnalysisContext& context){
    ClusterInputAdapterUtils::prepareSymmetryAwareClusterInputs(
        analysis,
        context,
        true,
        [](std::size_t, int structureType) {
            return structureType != LATTICE_OTHER;
        }
    );
}

}
