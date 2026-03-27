#include <volt/cli/common.h>
#include <volt/cna_service.h>
#include <volt/structures/crystal_structure_types.h>

#include <string>

using namespace Volt;
using namespace Volt::CLI;

namespace {

void showUsage(const std::string& name){
    printUsageHeader(name, "Volt - Common Neighbor Analysis");
    std::cerr
        << "  --crystalStructure <type>     Crystal structure. (FCC|BCC|HCP|CUBIC_DIAMOND|HEX_DIAMOND) [default: FCC]\n"
        << "  --dissolveSmallClusters       Mark small clusters as OTHER after building clusters.\n";
    printHelpOption();
}

}

int main(int argc, char* argv[]){
    if(argc < 2){
        showUsage(argv[0]);
        return 1;
    }

    std::string filename, outputBase;
    auto opts = parseArgs(argc, argv, filename, outputBase);

    if(hasOption(opts, "--help") || filename.empty()){
        showUsage(argv[0]);
        return filename.empty() ? 1 : 0;
    }

    initLogging("volt-common-neighbor-analysis");

    LammpsParser::Frame frame;
    if(!parseFrame(filename, frame)) return 1;

    outputBase = deriveOutputBase(filename, outputBase);
    spdlog::info("Output base: {}", outputBase);

    CommonNeighborAnalysisService analyzer;
    LatticeStructureType crystalStructure = LATTICE_FCC;
    const std::string crystalStructureOption = getString(opts, "--crystalStructure", "FCC");
    if(!parseLatticeStructureType(crystalStructureOption, crystalStructure)){
        spdlog::warn("Unknown crystal structure '{}', defaulting to FCC.", crystalStructureOption);
        crystalStructure = LATTICE_FCC;
    }
    analyzer.setInputCrystalStructure(crystalStructure);
    analyzer.setDissolveSmallClusters(hasOption(opts, "--dissolveSmallClusters"));

    spdlog::info("Starting common neighbor analysis...");
    json result = analyzer.compute(frame, outputBase);
    if(result.value("is_failed", false)){
        spdlog::error("Analysis failed: {}", result.value("error", "Unknown error"));
        return 1;
    }

    spdlog::info("Common neighbor analysis completed.");
    return 0;
}
