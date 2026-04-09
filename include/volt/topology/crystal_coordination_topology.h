#pragma once

#include <volt/core/volt.h>
#include <volt/topology/crystal_symmetry_topology.h>
#include <volt/topology/crystal_coordination_pattern.h>
#include <volt/analysis/nearest_neighbor_finder.h>
#include <volt/core/particle_property.h>
#include <volt/analysis/structure_analysis_context.h>
#include <volt/analysis/cna_local_structure.h>

namespace Volt{

class CoordinationStructures{
public:
    CoordinationStructures(
        ParticleProperty* structureTypes,
        LatticeStructureType inputCrystalType,
        bool identifyPlanarDefects,
        const SimulationCell& simCell
    );

    double determineLocalStructure(
        const NearestNeighborFinder& neighList,
        int particleIndex,
        int* outNeighborCount,
        int* outOrderedNeighborIndices = nullptr
    ) const;

    static void initializeStructures();

    void postProcessDiamondNeighbors(
        AnalysisContext& context,
        const NearestNeighborFinder& neighList
    ) const;

    static const LatticeStructureType getLatticeIdx(int s){
        switch(s){
            case StructureType::SC:
                return LATTICE_SC;
            case StructureType::FCC:
                return LATTICE_FCC;
            case StructureType::HCP:
                return LATTICE_HCP;
            case StructureType::BCC:
                return LATTICE_BCC;
            case StructureType::CUBIC_DIAMOND:
            case StructureType::CUBIC_DIAMOND_FIRST_NEIGH:
            case StructureType::CUBIC_DIAMOND_SECOND_NEIGH:
                return LATTICE_CUBIC_DIAMOND;
            case StructureType::HEX_DIAMOND:
            case StructureType::HEX_DIAMOND_FIRST_NEIGH:
            case StructureType::HEX_DIAMOND_SECOND_NEIGH:
                return LATTICE_HEX_DIAMOND;
            case StructureType::ICO:
            case StructureType::GRAPHENE:
                return LATTICE_OTHER;
            default:
                spdlog::warn("getLatticeIdx: unknown {}", s);
                return LATTICE_OTHER;
        }
    }

    static const CoordinationStructureType getCoordIdx(int s){
        switch(s){
            case StructureType::SC:
                return COORD_SC;
            case StructureType::FCC:
                return COORD_FCC;
            case StructureType::HCP:
                return COORD_HCP;
            case StructureType::BCC:
                return COORD_BCC;
            case StructureType::CUBIC_DIAMOND:
            case StructureType::CUBIC_DIAMOND_FIRST_NEIGH:
            case StructureType::CUBIC_DIAMOND_SECOND_NEIGH:
                return COORD_CUBIC_DIAMOND;
            case StructureType::HEX_DIAMOND:
            case StructureType::HEX_DIAMOND_FIRST_NEIGH:
            case StructureType::HEX_DIAMOND_SECOND_NEIGH:
                return COORD_HEX_DIAMOND;
            case StructureType::ICO:
            case StructureType::GRAPHENE:
                return COORD_OTHER;
            default:
                spdlog::warn("getCoordIdx: unknown {}", s);
                return COORD_OTHER;
        }
    }

    static const CoordinationStructure& getCoordStruct(int structureType){
        const auto index = static_cast<std::size_t>(CoordinationStructures::getCoordIdx(structureType));
        if(index >= _coordinationStructures.size()){
            return _coordinationStructures[static_cast<std::size_t>(COORD_OTHER)];
        }
        return _coordinationStructures[index];
    }

    static const LatticeStructure& getLatticeStruct(int structureType){
        const auto index = static_cast<std::size_t>(CoordinationStructures::getLatticeIdx(structureType));
        if(index >= _latticeStructures.size()){
            return _latticeStructures[static_cast<std::size_t>(LATTICE_OTHER)];
        }
        return _latticeStructures[index];
    }

    static const LatticeStructure& getLatticeStructByLatticeType(int latticeType){
        const auto index = static_cast<std::size_t>(latticeType);
        if(index >= _latticeStructures.size()){
            return _latticeStructures[static_cast<std::size_t>(LATTICE_OTHER)];
        }
        return _latticeStructures[index];
    }

    static std::vector<CoordinationStructure> _coordinationStructures;
    static std::vector<LatticeStructure> _latticeStructures;
    const SimulationCell& cell() const{
        return _simCell;
    }

    static void findNonCoplanarVectors(const CoordinationStructure& coordStruct, int nindices[3], Matrix3& tm1);

private:
    static void initializeDiamondStructure(
        int coordType,
        int latticeType,
        const Vector3* vectors,
        int numNeighbors,
        int totalVectors
    );

    static void initializeLatticeStructure(
        int latticeType,
        const Vector3* vectors,
        int totalVectors,
        CoordinationStructure* coordStruct
    );

    template <typename BondPredicate, typename SignatureFunction>
    static void initializeCoordinationStructure(
        int coordType,
        const Vector3* vectors,
        int numNeighbors,
        BondPredicate bondPredicate,
        SignatureFunction signatureFunction
    );

    static void initializeFCC();
    static void initializeSC();
    static void initializeHCP();
    static void initializeBCC();
    static void initializeCubicDiamond();
    static void initializeHexagonalDiamond();
    static void initializeOther();

    static void calculateSymmetryProducts(LatticeStructure& latticeStruct);
    static void generateSymmetryPermutations(LatticeStructure& latticeStruct);
    static void initializeSymmetryInformation();
    static void findCommonNeighborsForBond(CoordinationStructure& coordStruct, int neighborIndex);
    static void initializeCommonNeighbors();

    const SimulationCell& _simCell;
    ParticleProperty* _structureTypes;
    LatticeStructureType _inputCrystalType;
    bool _identifyPlanarDefects;
};

}
