#pragma once

#include <volt/analysis/analysis_context.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/analysis/cluster_builder.h>

#include <deque>

namespace Volt{

class CNAClusterBuilder final : public ClusterBuilder{
public:
    CNAClusterBuilder(
        StructureAnalysis& sa,
        AnalysisContext& context
    );

    void applyPreferredOrientation(Cluster* cluster);
    void build(bool dissolveSmallClusters = false);
    void grow(
        Cluster* cluster,
        std::deque<int>& atomsToVisit,
        Matrix_3<double>& orientationV,
        Matrix_3<double>& orientationW,
        int structureType
    );

private:
    void reorientAtomsToAlign();
    void connectClusters();
    void formSuperClusters();
    void processDefectCluster(Cluster* defectCluster);
    bool calculateMisorientation(
        int atomIndex,
        int neighbor,
        int neighborIndex,
        Matrix3& outTransition
    );
};

}
