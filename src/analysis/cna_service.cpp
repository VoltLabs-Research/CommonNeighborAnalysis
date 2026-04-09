#include <volt/analysis/cna_service.h>

#include <volt/analysis/reconstructed_analysis_pipeline.h>
#include <volt/analysis/structure_analysis_context.h>
#include <volt/analysis/cluster_graph_builder.h>
#include <volt/analysis/cluster_graph_io.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/analysis/cna_cluster_input_adapter.h>
#include <volt/analysis/cna_structure_analysis.h>
#include <volt/core/analysis_result.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/particle_property.h>
#include <volt/structures/crystal_structure_types.h>
#include <volt/utilities/json_utils.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <map>
#include <utility>

namespace Volt{

using namespace Volt::Particles;

namespace CnaServiceDetail{

std::string structureTypeNameForExport(int structureType){
    return structureTypeName(structureType);
}

json buildMainListing(const std::vector<int>& structureTypes){
    std::map<int, int> counts;

    for(int rawType : structureTypes){
        counts[rawType]++;
    }

    json listing = json::array();
    for(const auto& [structureType, count] : counts){
        if(count <= 0){
            continue;
        }

        listing.push_back({
            {"structure_type", structureType},
            {"structure_name", structureTypeNameForExport(structureType)},
            {"count", count}
        });
    }

    std::sort(listing.begin(), listing.end(), [](const json& lhs, const json& rhs){
        return lhs.value("structure_name", "") < rhs.value("structure_name", "");
    });

    return listing;
}

json buildPerAtomProperties(
    const LammpsParser::Frame& frame,
    const std::vector<int>& structureTypes,
    const ParticleProperty* clusterIds
){
    json perAtom = json::array();

    for(size_t atomIndex = 0; atomIndex < static_cast<size_t>(frame.natoms); ++atomIndex){
        const int structureType = structureTypes[atomIndex];

        json atom;
        atom["id"] = atomIndex < frame.ids.size()
            ? frame.ids[atomIndex]
            : static_cast<int>(atomIndex);
        atom["structure_type"] = structureType;
        atom["structure_name"] = structureTypeNameForExport(structureType);
        atom["cluster_id"] = clusterIds ? clusterIds->getInt(atomIndex) : 0;

        if(atomIndex < frame.positions.size()){
            const auto& pos = frame.positions[atomIndex];
            atom["pos"] = {pos.x(), pos.y(), pos.z()};
        }else{
            atom["pos"] = {0.0, 0.0, 0.0};
        }

        perAtom.push_back(std::move(atom));
    }

    return perAtom;
}

json buildAtomsExport(
    const LammpsParser::Frame& frame,
    const std::vector<int>& structureTypes,
    const ParticleProperty* clusterIds
){
    std::map<int, std::vector<size_t>> structureAtomIndices;
    for(size_t i = 0; i < static_cast<size_t>(frame.natoms); ++i){
        structureAtomIndices[structureTypes[i]].push_back(i);
    }

    std::vector<int> structureOrder;
    structureOrder.reserve(structureAtomIndices.size());
    for(const auto& [structureType, _] : structureAtomIndices){
        structureOrder.push_back(structureType);
    }
    std::sort(structureOrder.begin(), structureOrder.end(), [&](int a, int b){
        return structureTypeNameForExport(a) < structureTypeNameForExport(b);
    });

    json atomsByStructure = json::object();
    for(int st : structureOrder){
        json atomsArray = json::array();
        for(size_t atomIndex : structureAtomIndices[st]){
            const bool hasPosition = atomIndex < frame.positions.size();
            const auto atomId = atomIndex < frame.ids.size()
                ? frame.ids[atomIndex]
                : static_cast<int>(atomIndex);

            if(hasPosition){
                const auto& pos = frame.positions[atomIndex];
                atomsArray.push_back({
                    {"id", atomId},
                    {"cluster_id", clusterIds ? clusterIds->getInt(atomIndex) : 0},
                    {"pos", {pos.x(), pos.y(), pos.z()}}
                });
            }else{
                atomsArray.push_back({
                    {"id", atomId},
                    {"cluster_id", clusterIds ? clusterIds->getInt(atomIndex) : 0},
                    {"pos", {0.0, 0.0, 0.0}}
                });
            }
        }
        atomsByStructure[structureTypeNameForExport(st)] = atomsArray;
    }

    json exportWrapper;
    exportWrapper["export"] = json::object();
    exportWrapper["export"]["AtomisticExporter"] = atomsByStructure;
    return exportWrapper;
}

}

CommonNeighborAnalysisService::CommonNeighborAnalysisService()
    : _inputCrystalStructure(LATTICE_FCC)
    , _dissolveSmallClusters(false){}

void CommonNeighborAnalysisService::setInputCrystalStructure(LatticeStructureType structureType){
    _inputCrystalStructure = structureType;
}

void CommonNeighborAnalysisService::setDissolveSmallClusters(bool dissolveSmallClusters){
    _dissolveSmallClusters = dissolveSmallClusters;
}

json CommonNeighborAnalysisService::compute(
    const LammpsParser::Frame& frame,
    const std::string& outputBase,
    const std::string& inputDumpPath
){
    const std::string annotatedDumpPath = outputBase.empty()
        ? inputDumpPath
        : outputBase + "_annotated.dump";

    if(_inputCrystalStructure == LATTICE_SC){
        return AnalysisResult::failure("CNA does not support SC. Use PTM for this crystal.");
    }

    std::string frameError;
    auto session = AnalysisPipelineUtils::prepareAnalysisSession(
        frame,
        _inputCrystalStructure,
        &frameError
    );
    if(!session){
        return AnalysisResult::failure(frameError);
    }
    AnalysisContext& context = session->context;

    try{
        StructureAnalysis analysis(context);
        identifyStructuresCNA(analysis);
        CNAClusterInputAdapter clusterInputAdapter;
        clusterInputAdapter.prepare(analysis, context);
        ClusterBuilder clusterBuilder(analysis, context);
        clusterBuilder.build(_dissolveSmallClusters);

        std::vector<int> atomStructureTypes(
            static_cast<size_t>(frame.natoms),
            static_cast<int>(StructureType::OTHER)
        );
        for(int atomIndex = 0; atomIndex < frame.natoms; ++atomIndex){
            atomStructureTypes[static_cast<size_t>(atomIndex)] = context.structureTypes->getInt(atomIndex);
        }

        json result;
        result["main_listing"] = CnaServiceDetail::buildMainListing(atomStructureTypes);
        result["per-atom-properties"] = CnaServiceDetail::buildPerAtomProperties(
            frame,
            atomStructureTypes,
            context.atomClusters.get()
        );

        if(!outputBase.empty()){
            const std::string msgpackPath = outputBase + "_cna_analysis.msgpack";
            if(!JsonUtils::writeJsonMsgpackToFile(result, msgpackPath, false)){
                return AnalysisResult::failure("Failed to write " + msgpackPath);
            }

            const std::string atomsPath = outputBase + "_atoms.msgpack";
            if(!JsonUtils::writeJsonMsgpackToFile(
                CnaServiceDetail::buildAtomsExport(frame, atomStructureTypes, context.atomClusters.get()),
                atomsPath,
                false
            )){
                return AnalysisResult::failure("Failed to write " + atomsPath);
            }
        }

        if(!AnalysisPipelineUtils::appendClusterOutputs(
            frame,
            outputBase,
            annotatedDumpPath,
            context,
            analysis,
            result,
            &frameError
        )){
            return AnalysisResult::failure(frameError);
        }

        result["is_failed"] = false;
        return result;
    }catch(const std::exception& error){
        return AnalysisResult::failure(std::string("CNA analysis failed: ") + error.what());
    }
}

}
