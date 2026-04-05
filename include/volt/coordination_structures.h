#pragma once

#include <volt/core/volt.h>
#include <volt/lattice_structure.h>
#include <volt/coordination_structure.h>
#include <volt/analysis/nearest_neighbor_finder.h>
#include <volt/core/particle_property.h>
#include <volt/structures/lattice_vectors.h>
#include <volt/analysis/analysis_context.h>
#include <volt/cna_local_structure_utils.h>

namespace Volt{
    
class CoordinationStructures{
public:
    CoordinationStructures(ParticleProperty* structureTypes, LatticeStructureType inputCrystalType, bool identifyPlanarDefects, const SimulationCell& simCell);

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
            case StructureType::SC:  return LATTICE_SC;
            case StructureType::FCC: return LATTICE_FCC;
            case StructureType::HCP: return LATTICE_HCP;
            case StructureType::BCC: return LATTICE_BCC;

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
            case StructureType::SC:  return COORD_SC;
            case StructureType::FCC: return COORD_FCC;
            case StructureType::HCP: return COORD_HCP;
            case StructureType::BCC: return COORD_BCC;

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
        return _coordinationStructures[CoordinationStructures::getCoordIdx(structureType)];
    }

    static const LatticeStructure& getLatticeStruct(int structureType){
        return _latticeStructures[CoordinationStructures::getLatticeIdx(structureType)];
    }

    static CoordinationStructure _coordinationStructures[NUM_COORD_TYPES];
    static LatticeStructure _latticeStructures[NUM_LATTICE_TYPES];
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
