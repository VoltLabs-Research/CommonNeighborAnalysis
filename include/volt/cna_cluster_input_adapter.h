#pragma once

#include <volt/analysis/cluster_input_adapter.h>

namespace Volt{

class CNAClusterInputAdapter final : public ClusterInputAdapter{
public:
    void prepare(StructureAnalysis& analysis, AnalysisContext& context) override;
};

}
