#pragma once

#include <volt/core/volt.h>
#include <nlohmann/json.hpp>
#include <volt/core/lammps_parser.h>
#include <volt/structures/crystal_structure_types.h>
#include <string>

namespace Volt{

using json = nlohmann::json;

class CommonNeighborAnalysisService{
public:
    CommonNeighborAnalysisService();

    void setInputCrystalStructure(LatticeStructureType structureType);

    json compute(
        const LammpsParser::Frame& frame,
        const std::string& outputBase = ""
    );

private:
    LatticeStructureType _inputCrystalStructure;
};

}
