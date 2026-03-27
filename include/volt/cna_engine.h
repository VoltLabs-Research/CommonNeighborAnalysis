#pragma once

#include <volt/core/volt.h>
#include <volt/analysis/analysis_context.h>
#include <volt/cna.h>
#include <volt/analysis/nearest_neighbor_finder.h>
#include <volt/coordination_structures.h>
#include <volt/core/lammps_parser.h>
#include <nlohmann/json.hpp>
#include <array>
#include <map>

namespace Volt{

using json = nlohmann::json;

class CommonNeighborAnalysisEngine{
public:
    CommonNeighborAnalysisEngine(AnalysisContext& context, bool identifyPlanarDefects);

    void perform();

    json buildMainListing() const;
    json getPerAtomProperties(const LammpsParser::Frame& frame) const;

    std::string getStructureTypeName(int structureType) const;

private:
    int coordinationNumber() const;
    double computeLocalCutoff(
        const NearestNeighborFinder& neighList,
        const NearestNeighborFinder::Query<MAX_NEIGHBORS>& neighQuery,
        int numNeighbors,
        int coordinationNumber,
        int particleIndex,
        int* neighborIndices,
        Vector3* neighborVectors,
        NeighborBondArray& neighborArray
    ) const;

    double determineLocalStructure(
        const NearestNeighborFinder& neighList,
        int particleIndex
    ) const;

    void storeNeighborOrdering(
        const std::array<int, MAX_NEIGHBORS>& neighborMapping,
        int coordinationCount,
        StructureType atomStructure,
        int particleIndex
    ) const;

    AnalysisContext& _context;
};

}
