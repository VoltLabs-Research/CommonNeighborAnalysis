#include <volt/cli/common.h>
#include <volt/cna_service.h>
#include <volt/structures/crystal_structure_types.h>

#include <string>

using namespace Volt;
using namespace Volt::CLI;

namespace {

LatticeStructureType parseCrystalStructure(const std::string& value){
    if(value == "FCC") return LATTICE_FCC;
    if(value == "BCC") return LATTICE_BCC;
    if(value == "HCP") return LATTICE_HCP;
    if(value == "SC") return LATTICE_SC;
    if(value == "CUBIC_DIAMOND") return LATTICE_CUBIC_DIAMOND;
    if(value == "HEX_DIAMOND") return LATTICE_HEX_DIAMOND;
    spdlog::warn("Unknown crystal structure '{}', defaulting to FCC.", value);
    return LATTICE_FCC;
}

void showUsage(const std::string& name){
    printUsageHeader(name, "Volt - Common Neighbor Analysis");
    std::cerr
        << "  --crystalStructure <type>     Crystal structure. (FCC|BCC|HCP|SC|CUBIC_DIAMOND|HEX_DIAMOND) [default: FCC]\n"
        << "  --threads <int>               Max worker threads (TBB/OMP). [default: auto]\n";
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

    auto parallel = initParallelism(opts, false);
    initLogging("volt-common-neighbor-analysis", parallel.threads);

    LammpsParser::Frame frame;
    if(!parseFrame(filename, frame)) return 1;

    outputBase = deriveOutputBase(filename, outputBase);
    spdlog::info("Output base: {}", outputBase);

    CommonNeighborAnalysisService analyzer;
    analyzer.setInputCrystalStructure(parseCrystalStructure(getString(opts, "--crystalStructure", "FCC")));

    spdlog::info("Starting common neighbor analysis...");
    json result = analyzer.compute(frame, outputBase);
    if(result.value("is_failed", false)){
        spdlog::error("Analysis failed: {}", result.value("error", "Unknown error"));
        return 1;
    }

    spdlog::info("Common neighbor analysis completed.");
    return 0;
}
