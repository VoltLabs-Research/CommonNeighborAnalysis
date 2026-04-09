#pragma once

#include <volt/core/volt.h>
#include <volt/topology/crystal_symmetry_topology.h>
#include <volt/topology/crystal_coordination_pattern.h>
#include <volt/analysis/nearest_neighbor_finder.h>
#include <volt/core/particle_property.h>
#include <volt/structures/crystal_topology_registry.h>
#include <volt/analysis/structure_analysis_context.h>
#include <volt/analysis/cna_local_structure.h>

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
        const auto* topology = crystalTopologyByStructureType(s);
        if(!topology || topology->latticeType <= 0){
            if(s != StructureType::OTHER){
                spdlog::warn("getLatticeIdx: unknown {}", s);
            }
            return LATTICE_OTHER;
        }
        return static_cast<LatticeStructureType>(topology->latticeType);
    }

    static const CoordinationStructureType getCoordIdx(int s){
        const auto* topology = crystalTopologyByStructureType(s);
        if(!topology || topology->coordinationType <= 0){
            if(s != StructureType::OTHER){
                spdlog::warn("getCoordIdx: unknown {}", s);
            }
            return COORD_OTHER;
        }
        return static_cast<CoordinationStructureType>(topology->coordinationType);
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
    static void initializeOther();
    static void initializeFromRegistry();
    
    static void calculateSymmetryProducts(LatticeStructure& latticeStruct);
    static void initializeSymmetryInformation();
        
    const SimulationCell& _simCell;
    ParticleProperty* _structureTypes;
    LatticeStructureType _inputCrystalType;
    bool _identifyPlanarDefects;
};

}
