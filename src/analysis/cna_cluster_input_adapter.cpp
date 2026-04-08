#include <volt/analysis/cna_cluster_input_adapter.h>
#include <volt/analysis/cluster_input_preparation.h>

namespace Volt{

void CNAClusterInputAdapter::prepare(StructureAnalysis& analysis, AnalysisContext& context){
    ClusterInputAdapterUtils::prepareSymmetryAwareClusterInputs(
        analysis,
        context,
        true,
        [&](std::size_t atomIndex, int structureType) {
            if(structureType == LATTICE_OTHER){
                return false;
            }
            if(analysis.numberOfNeighbors(static_cast<int>(atomIndex)) == 0){
                return false;
            }
            return true;
        }
    );
}

}
