#include <volt/cna_service.h>

#include <volt/analysis/analysis_context.h>
#include <volt/analysis/cluster_graph_export.h>
#include <volt/analysis/structure_analysis.h>
#include <volt/cna_cluster_builder.h>
#include <volt/core/analysis_result.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/particle_property.h>
#include <volt/structures/crystal_structure_types.h>
#include <volt/utilities/json_utils.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <utility>

namespace Volt{

using namespace Volt::Particles;

namespace{

std::string structureTypeNameForExport(int structureType){
    switch(static_cast<StructureType>(structureType)){
        case StructureType::SC:
            return "SC";
        case StructureType::FCC:
            return "FCC";
        case StructureType::HCP:
            return "HCP";
        case StructureType::BCC:
            return "BCC";
        case StructureType::CUBIC_DIAMOND:
            return "CUBIC_DIAMOND";
        case StructureType::HEX_DIAMOND:
            return "HEX_DIAMOND";
        case StructureType::ICO:
            return "ICO";
        case StructureType::GRAPHENE:
            return "GRAPHENE";
        case StructureType::CUBIC_DIAMOND_FIRST_NEIGH:
            return "CUBIC_DIAMOND_FIRST_NEIGH";
        case StructureType::CUBIC_DIAMOND_SECOND_NEIGH:
            return "CUBIC_DIAMOND_SECOND_NEIGH";
        case StructureType::HEX_DIAMOND_FIRST_NEIGH:
            return "HEX_DIAMOND_FIRST_NEIGH";
        case StructureType::HEX_DIAMOND_SECOND_NEIGH:
            return "HEX_DIAMOND_SECOND_NEIGH";
        case StructureType::OTHER:
        case StructureType::NUM_STRUCTURE_TYPES:
        default:
            return "OTHER";
    }
}

json buildMainListing(const std::vector<int>& structureTypes){
    constexpr int K = static_cast<int>(StructureType::NUM_STRUCTURE_TYPES);
    std::vector<int> counts(static_cast<size_t>(K), 0);

    for(int rawType : structureTypes){
        const int safeType = (0 <= rawType && rawType < K)
            ? rawType
            : static_cast<int>(StructureType::OTHER);
        counts[static_cast<size_t>(safeType)]++;
    }

    json listing = json::array();
    for(int structureType = 0; structureType < K; ++structureType){
        const int count = counts[static_cast<size_t>(structureType)];
        if(count == 0){
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
        const int rawType = structureTypes[atomIndex];
        const int structureType = (0 <= rawType && rawType < static_cast<int>(StructureType::NUM_STRUCTURE_TYPES))
            ? rawType
            : static_cast<int>(StructureType::OTHER);

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
    constexpr int K = static_cast<int>(StructureType::NUM_STRUCTURE_TYPES);

    std::vector<std::string> names(K);
    for(int st = 0; st < K; ++st){
        names[st] = structureTypeNameForExport(st);
    }

    std::vector<std::vector<size_t>> structureAtomIndices(K);
    for(size_t i = 0; i < static_cast<size_t>(frame.natoms); ++i){
        const int raw = structureTypes[i];
        const int st = (0 <= raw && raw < K) ? raw : static_cast<int>(StructureType::OTHER);
        structureAtomIndices[static_cast<size_t>(st)].push_back(i);
    }

    std::vector<int> structureOrder;
    structureOrder.reserve(K);
    for(int st = 0; st < K; ++st){
        if(!structureAtomIndices[static_cast<size_t>(st)].empty()){
            structureOrder.push_back(st);
        }
    }

    std::sort(structureOrder.begin(), structureOrder.end(), [&](int a, int b){
        return names[a] < names[b];
    });

    json atomsByStructure = json::object();
    for(int st : structureOrder){
        json atomsArray = json::array();
        for(size_t atomIndex : structureAtomIndices[static_cast<size_t>(st)]){
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
        atomsByStructure[names[st]] = atomsArray;
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
    const std::string& outputBase
){
    if(_inputCrystalStructure == LATTICE_SC){
        return AnalysisResult::failure("CNA does not support SC. Use PTM for simple cubic crystals.");
    }

    FrameAdapter::PreparedAnalysisInput prepared;
    std::string frameError;
    if(!FrameAdapter::prepareAnalysisInput(frame, prepared, &frameError)){
        return AnalysisResult::failure(frameError);
    }

    auto positions = std::move(prepared.positions);

    auto structureTypes = std::make_unique<ParticleProperty>(frame.natoms, DataType::Int, 1, 0, true);
    std::vector<Matrix3> preferredOrientations{Matrix3::Identity()};

    AnalysisContext context(
        positions.get(),
        frame.simulationCell,
        _inputCrystalStructure,
        nullptr,
        structureTypes.get(),
        std::move(preferredOrientations)
    );

    try{
        StructureAnalysis analysis(
            context
        );
        analysis.identifyStructuresCNA();
        CNAClusterBuilder clusterBuilder(analysis, context);
        clusterBuilder.build(_dissolveSmallClusters);

        std::vector<int> atomStructureTypes(
            static_cast<size_t>(frame.natoms),
            static_cast<int>(StructureType::OTHER)
        );
        for(int atomIndex = 0; atomIndex < frame.natoms; ++atomIndex){
            atomStructureTypes[static_cast<size_t>(atomIndex)] = context.structureTypes->getInt(atomIndex);
        }

        json result;
        result["main_listing"] = buildMainListing(atomStructureTypes);
        result["per-atom-properties"] = buildPerAtomProperties(
            frame,
            atomStructureTypes,
            context.atomClusters.get()
        );

        if(!outputBase.empty()){
            const std::string dumpPath = outputBase + "_annotated.dump";
            if(!context.writeDumpWithContext(frame, dumpPath, &analysis)){
                return AnalysisResult::failure("Failed to write " + dumpPath);
            }

            const std::string msgpackPath = outputBase + "_common_neighbor_analysis.msgpack";
            if(!JsonUtils::writeJsonMsgpackToFile(result, msgpackPath, false)){
                return AnalysisResult::failure("Failed to write " + msgpackPath);
            }

            const std::string atomsPath = outputBase + "_atoms.msgpack";
            if(!JsonUtils::writeJsonMsgpackToFile(
                buildAtomsExport(frame, atomStructureTypes, context.atomClusters.get()),
                atomsPath,
                false
            )){
                return AnalysisResult::failure("Failed to write " + atomsPath);
            }

            ClusterGraphExportPaths clusterGraphPaths;
            if(!exportClusterGraph(analysis.clusterGraph(), outputBase, &clusterGraphPaths)){
                return AnalysisResult::failure("Failed to export cluster graph tables");
            }

            result["annotated_dump"] = dumpPath;
            result["clusters_table"] = clusterGraphPaths.clustersTablePath;
            result["cluster_transitions_table"] = clusterGraphPaths.clusterTransitionsTablePath;
        }

        result["is_failed"] = false;
        return result;
    }catch(const std::exception& error){
        return AnalysisResult::failure(std::string("CNA analysis failed: ") + error.what());
    }
}

}
