#include <volt/topology/crystal_coordination_topology.h>
#include <volt/analysis/structure_analysis_context.h>
#include <volt/analysis/cna_classifier.h>
#include <volt/analysis/crystal_symmetry_utils.h>

#include <stdexcept>

namespace Volt{

// Contains the known coordination structures.
std::vector<CoordinationStructure> CoordinationStructures::_coordinationStructures;

// Contains the known lattice types.
std::vector<LatticeStructure> CoordinationStructures::_latticeStructures;

CoordinationStructures::CoordinationStructures(
	ParticleProperty* structureTypes,
	LatticeStructureType inputCrystalType,
	bool identifyPlanarDefects,
	const SimulationCell& simCell
	) : _structureTypes(structureTypes)
		, _inputCrystalType(inputCrystalType)
		, _identifyPlanarDefects(identifyPlanarDefects)
		, _simCell(simCell){}

// Determines the coordination structure of a particle
double CoordinationStructures::determineLocalStructure(
	const NearestNeighborFinder& neighList, 
	int particleIndex,
	int* outNeighborCount,
	int* outOrderedNeighborIndices
) const { 
    assert(_structureTypes->getInt(particleIndex) == COORD_OTHER);
    CnaLocalStructureUtils::LocalStructureMatch match;
    if(!CnaLocalStructureUtils::determineLocalStructure(
        _inputCrystalType,
        _identifyPlanarDefects,
        neighList,
        particleIndex,
        _coordinationStructures.data(),
        match
    )){
        return 0.0;
    }

    _structureTypes->setInt(particleIndex, static_cast<int>(match.structureType));

	if(outNeighborCount){
		*outNeighborCount = match.coordinationNumber;
	}

	if(outOrderedNeighborIndices){
		for(int i = 0; i < match.coordinationNumber; ++i){
			outOrderedNeighborIndices[i] = match.neighborIndices[match.neighborMapping[i]];
		}
	}

	return match.localCutoff;
}

void CoordinationStructures::postProcessDiamondNeighbors(
    AnalysisContext& context,
    const NearestNeighborFinder& neighList
) const{
    if(_inputCrystalType != LATTICE_CUBIC_DIAMOND && _inputCrystalType != LATTICE_HEX_DIAMOND) return;
    
    const size_t N = _structureTypes->size();
    std::vector<int> newStructureTypes(N);

    for(size_t i = 0; i < N; ++i){
        newStructureTypes[i] = _structureTypes->getInt(i);
    }
    
    NearestNeighborFinder firstNeighborFinder(4);
    firstNeighborFinder.prepare(context.positions, context.simCell, nullptr);
    
    // Mark first neighbors of diamond atoms
    for(size_t i = 0; i < N; ++i){
        int currentType = _structureTypes->getInt(i);
        if(currentType != static_cast<int>(StructureType::CUBIC_DIAMOND) && 
           currentType != static_cast<int>(StructureType::HEX_DIAMOND)) continue;
           
        NearestNeighborFinder::Query<4> query(firstNeighborFinder);
        // false = do not include self
        query.findNeighbors(i, false); 
        
        StructureType firstNeighType = (currentType == static_cast<int>(StructureType::CUBIC_DIAMOND)) 
            ? StructureType::CUBIC_DIAMOND_FIRST_NEIGH 
            : StructureType::HEX_DIAMOND_FIRST_NEIGH;
            
        for(const auto& neighbor : query.results()){
            if(_structureTypes->getInt(neighbor.index) == static_cast<int>(StructureType::OTHER)){
                newStructureTypes[neighbor.index] = static_cast<int>(firstNeighType);
            }
        }
    }
    
    // Apply changes for first neighbors
    for(size_t i = 0; i < N; ++i){
        _structureTypes->setInt(i, newStructureTypes[i]);
    }
    
    // Mark second neighbors from first neighbors
    for(size_t i = 0; i < N; ++i){
        int currentType = _structureTypes->getInt(i);
        if(currentType != static_cast<int>(StructureType::CUBIC_DIAMOND_FIRST_NEIGH) && 
           currentType != static_cast<int>(StructureType::HEX_DIAMOND_FIRST_NEIGH)) continue;
           
        NearestNeighborFinder::Query<4> query(firstNeighborFinder);
        query.findNeighbors(i, false);
        
        StructureType secondNeighType = (currentType == static_cast<int>(StructureType::CUBIC_DIAMOND_FIRST_NEIGH)) 
            ? StructureType::CUBIC_DIAMOND_SECOND_NEIGH 
            : StructureType::HEX_DIAMOND_SECOND_NEIGH;
            
        for(const auto& neighbor : query.results()){
            if(_structureTypes->getInt(neighbor.index) == static_cast<int>(StructureType::OTHER)){
                newStructureTypes[neighbor.index] = static_cast<int>(secondNeighType);
            }
        }
    }
    
    for(size_t i = 0; i < N; ++i){
        _structureTypes->setInt(i, newStructureTypes[i]);
    }
}

void populateCoordinationStructureFromTopology(
    CoordinationStructure& coordStruct,
    const CrystalTopologyEntry& entry
){
    coordStruct.numNeighbors = entry.coordinationNumber;
    coordStruct.latticeVectors.assign(
        entry.latticeVectors.begin(),
        entry.latticeVectors.begin() + entry.coordinationNumber
    );

    for(int row = 0; row < entry.coordinationNumber; ++row){
        coordStruct.neighborArray.neighborArray[static_cast<std::size_t>(row)] =
            entry.neighborBondRows[static_cast<std::size_t>(row)];
        coordStruct.cnaSignatures[row] = row < static_cast<int>(entry.cnaSignatureCodes.size())
            ? entry.cnaSignatureCodes[static_cast<std::size_t>(row)]
            : 0;
        coordStruct.commonNeighbors[row][0] =
            entry.commonNeighbors[static_cast<std::size_t>(row)][0];
        coordStruct.commonNeighbors[row][1] =
            entry.commonNeighbors[static_cast<std::size_t>(row)][1];
    }
}

void populateLatticeStructureFromTopology(
    LatticeStructure& latticeStruct,
    const CrystalTopologyEntry& entry,
    CoordinationStructure* coordStruct
){
    latticeStruct.coordStructure = coordStruct;
    latticeStruct.latticeVectors = entry.latticeVectors;
    latticeStruct.primitiveCell = entry.primitiveCell;
    latticeStruct.primitiveCellInverse = entry.primitiveCellInverse;
    latticeStruct.maxNeighbors = coordStruct->numNeighbors;
    latticeStruct.permutations.clear();
    latticeStruct.permutations.reserve(entry.symmetries.size());
    for(const auto& symmetry : entry.symmetries){
        SymmetryPermutation permutation;
        permutation.transformation = symmetry.transformation;
        permutation.permutation.fill(-1);
        for(std::size_t slot = 0; slot < symmetry.permutation.size() && slot < permutation.permutation.size(); ++slot){
            permutation.permutation[slot] = symmetry.permutation[slot];
        }
        latticeStruct.permutations.push_back(std::move(permutation));
    }
}

void CoordinationStructures::initializeFromRegistry(){
    std::size_t maxCoordinationType = static_cast<std::size_t>(COORD_OTHER);
    std::size_t maxLatticeType = static_cast<std::size_t>(LATTICE_OTHER);
    for(const CrystalTopologyEntry& entry : crystalTopologyRegistry().entries()){
        if(entry.coordinationType > 0){
            maxCoordinationType = std::max(maxCoordinationType, static_cast<std::size_t>(entry.coordinationType));
        }
        if(entry.latticeType > 0){
            maxLatticeType = std::max(maxLatticeType, static_cast<std::size_t>(entry.latticeType));
        }
    }

    if(_coordinationStructures.size() <= maxCoordinationType){
        _coordinationStructures.resize(maxCoordinationType + 1);
    }
    if(_latticeStructures.size() <= maxLatticeType){
        _latticeStructures.resize(maxLatticeType + 1);
    }

    for(const CrystalTopologyEntry& entry : crystalTopologyRegistry().entries()){
        if(entry.coordinationType <= 0 || entry.latticeType <= 0){
            continue;
        }
        populateCoordinationStructureFromTopology(
            _coordinationStructures[static_cast<std::size_t>(entry.coordinationType)],
            entry
        );
        populateLatticeStructureFromTopology(
            _latticeStructures[static_cast<std::size_t>(entry.latticeType)],
            entry,
            &_coordinationStructures[static_cast<std::size_t>(entry.coordinationType)]
        );
    }
}

void CoordinationStructures::initializeOther(){
    if(_coordinationStructures.empty()){
        _coordinationStructures.resize(static_cast<std::size_t>(COORD_OTHER) + 1);
    }
    if(_latticeStructures.empty()){
        _latticeStructures.resize(static_cast<std::size_t>(LATTICE_OTHER) + 1);
    }

    _coordinationStructures[static_cast<std::size_t>(COORD_OTHER)].numNeighbors = 0;
    _latticeStructures[static_cast<std::size_t>(LATTICE_OTHER)].coordStructure =
        &_coordinationStructures[static_cast<std::size_t>(COORD_OTHER)];
    _latticeStructures[static_cast<std::size_t>(LATTICE_OTHER)].primitiveCell.setZero();
    _latticeStructures[static_cast<std::size_t>(LATTICE_OTHER)].primitiveCellInverse.setZero();
    _latticeStructures[static_cast<std::size_t>(LATTICE_OTHER)].maxNeighbors = 0;
    _latticeStructures[static_cast<std::size_t>(LATTICE_OTHER)].latticeVectors.clear();
    _latticeStructures[static_cast<std::size_t>(LATTICE_OTHER)].permutations.clear();
}

void CoordinationStructures::initializeSymmetryInformation(){
    for(auto latticeStruct = std::begin(_latticeStructures);
        latticeStruct != std::end(_latticeStructures); ++latticeStruct){
        
        if(latticeStruct->latticeVectors.empty()) continue;
        
        latticeStruct->primitiveCellInverse = latticeStruct->primitiveCell.inverse();
        if(latticeStruct->permutations.empty()){
            throw std::runtime_error("Missing explicit symmetry_permutations in topology metadata.");
        }
        calculateSymmetryProducts(*latticeStruct);
    }
}

void CoordinationStructures::findNonCoplanarVectors(const CoordinationStructure& coordStruct, int nindices[3], Matrix3& tm1){
    AnalysisSymmetryUtils::findNonCoplanarVectors(coordStruct.latticeVectors, coordStruct.numNeighbors, nindices, tm1);
}

void CoordinationStructures::calculateSymmetryProducts(LatticeStructure& latticeStruct){
    AnalysisSymmetryUtils::calculateSymmetryProducts(latticeStruct.permutations);
}

void CoordinationStructures::initializeStructures(){
	initializeOther();
    initializeFromRegistry();
	initializeSymmetryInformation();
}

}
