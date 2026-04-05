#include <volt/topology/crystal_coordination_topology.h>
#include <volt/analysis/structure_analysis_context.h>
#include <volt/analysis/cna_classifier.h>
#include <volt/analysis/crystal_symmetry_utils.h>

namespace Volt{

// Contains the known coordination structures.
CoordinationStructure CoordinationStructures::_coordinationStructures[NUM_COORD_TYPES];

// Contains the known lattice types.
LatticeStructure CoordinationStructures::_latticeStructures[NUM_LATTICE_TYPES];

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
        _coordinationStructures,
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

void CoordinationStructures::initializeFCC(){
	initializeCoordinationStructure(COORD_FCC, FCC_VECTORS, 12, [&](const Vector3& v1, const Vector3& v2){
		return (v1 - v2).length() < (sqrt(0.5f) + 1.0) * 0.5;
	}, [](int ni){ return 0; });

	initializeLatticeStructure(LATTICE_FCC, FCC_VECTORS, 12, &_coordinationStructures[COORD_FCC]);
    _latticeStructures[LATTICE_FCC].primitiveCell.column(0) = FCC_PRIMITIVE_CELL[0];
    _latticeStructures[LATTICE_FCC].primitiveCell.column(1) = FCC_PRIMITIVE_CELL[1];
    _latticeStructures[LATTICE_FCC].primitiveCell.column(2) = FCC_PRIMITIVE_CELL[2];
}

void CoordinationStructures::initializeHCP(){
	initializeCoordinationStructure(COORD_HCP, HCP_VECTORS, 12, [&](const Vector3& v1, const Vector3& v2){
		return (v1 - v2).length() < (sqrt(0.5) + 1.0) * 0.5;
	}, [&](int ni){ return (HCP_VECTORS[ni].z() == 0) ? 1 : 0; });

	initializeLatticeStructure(LATTICE_HCP, HCP_VECTORS, 18, &_coordinationStructures[COORD_HCP]);
    _latticeStructures[LATTICE_HCP].primitiveCell.column(0) = HCP_PRIMITIVE_CELL[0];
    _latticeStructures[LATTICE_HCP].primitiveCell.column(1) = HCP_PRIMITIVE_CELL[1];
    _latticeStructures[LATTICE_HCP].primitiveCell.column(2) = HCP_PRIMITIVE_CELL[2];
}

void CoordinationStructures::initializeBCC(){
	initializeCoordinationStructure(COORD_BCC, BCC_VECTORS, 14, [&](const Vector3& v1, const Vector3& v2){
		return (v1 - v2).length() < (double(1) + sqrt(double(2))) * double(0.5);
	}, [](int ni) { return (ni < 8) ? 0 : 1; });

	initializeLatticeStructure(LATTICE_BCC, BCC_VECTORS, 14, &_coordinationStructures[COORD_BCC]);
    _latticeStructures[LATTICE_BCC].primitiveCell.column(0) = BCC_PRIMITIVE_CELL[0];
    _latticeStructures[LATTICE_BCC].primitiveCell.column(1) = BCC_PRIMITIVE_CELL[1];
    _latticeStructures[LATTICE_BCC].primitiveCell.column(2) = BCC_PRIMITIVE_CELL[2];
}

void CoordinationStructures::initializeCubicDiamond(){
    initializeDiamondStructure(COORD_CUBIC_DIAMOND, LATTICE_CUBIC_DIAMOND, DIAMOND_CUBIC_VECTORS, 16, 20);
    
    _latticeStructures[LATTICE_CUBIC_DIAMOND].primitiveCell.column(0) = CUBIC_DIAMOND_PRIMITIVE_CELL[0];
    _latticeStructures[LATTICE_CUBIC_DIAMOND].primitiveCell.column(1) = CUBIC_DIAMOND_PRIMITIVE_CELL[1];
    _latticeStructures[LATTICE_CUBIC_DIAMOND].primitiveCell.column(2) = CUBIC_DIAMOND_PRIMITIVE_CELL[2];
}

void CoordinationStructures::initializeHexagonalDiamond(){
	initializeDiamondStructure(COORD_HEX_DIAMOND, LATTICE_HEX_DIAMOND, DIAMOND_HEX_VECTORS, 16, 32);
    
    _latticeStructures[LATTICE_HEX_DIAMOND].primitiveCell.column(0) = HEXAGONAL_DIAMOND_PRIMITIVE_CELL[0];
    _latticeStructures[LATTICE_HEX_DIAMOND].primitiveCell.column(1) = HEXAGONAL_DIAMOND_PRIMITIVE_CELL[1];
    _latticeStructures[LATTICE_HEX_DIAMOND].primitiveCell.column(2) = HEXAGONAL_DIAMOND_PRIMITIVE_CELL[2];
}

void CoordinationStructures::initializeSC(){
    initializeCoordinationStructure(COORD_SC, SC_VECTORS, 6, [](const Vector3& v1, const Vector3& v2){
        const double dot = v1.dot(v2);
        // 90°
        if(std::abs(dot) < EPSILON) return true;
        // opposites
        if((v1 + v2).squaredLength() < EPSILON) return true;
        return false;
    }, [](int /* ni */){ return 0; });

    initializeLatticeStructure(LATTICE_SC, SC_VECTORS, 6, &_coordinationStructures[COORD_SC]);
    _latticeStructures[LATTICE_SC].primitiveCell.column(0) = SC_PRIMITIVE_CELL[0];
    _latticeStructures[LATTICE_SC].primitiveCell.column(1) = SC_PRIMITIVE_CELL[1];
    _latticeStructures[LATTICE_SC].primitiveCell.column(2) = SC_PRIMITIVE_CELL[2];

    // permutations in all 6 directions
    _latticeStructures[LATTICE_SC].permutations.clear();
    for(const auto& R : AnalysisSymmetryUtils::cubicSymmetryRotations()){
        SymmetryPermutation sp;
        sp.transformation = R;
        for(int i = 0; i < 6; ++i){
            const Vector3 vec = R * SC_VECTORS[i];
            // match with index j such that R * vi == vj
            for(int j = 0; j < 6; ++j){
                if(vec.equals(SC_VECTORS[j])){
                    sp.permutation[i] = j;
                    break;
                }
            }
        }

        _latticeStructures[LATTICE_SC].permutations.push_back(sp);
    }
}   

void CoordinationStructures::initializeOther(){
    _coordinationStructures[COORD_OTHER].numNeighbors = 0;
    _latticeStructures[LATTICE_OTHER].coordStructure = &_coordinationStructures[COORD_OTHER];
    _latticeStructures[LATTICE_OTHER].primitiveCell.setZero();
    _latticeStructures[LATTICE_OTHER].primitiveCellInverse.setZero();
    _latticeStructures[LATTICE_OTHER].maxNeighbors = 0;
}

template <typename BondPredicate, typename SignatureFunction>
void CoordinationStructures::initializeCoordinationStructure(
	int coordType,
	const Vector3* vectors,
	int numNeighbors,
	BondPredicate bondPred,
	SignatureFunction sigFunc
){
	_coordinationStructures[coordType].numNeighbors = numNeighbors;
	for(int ni1 = 0; ni1 < numNeighbors; ni1++){
        _coordinationStructures[coordType].neighborArray.setNeighborBond(ni1, ni1, false);
        for(int ni2 = ni1 + 1; ni2 < numNeighbors; ni2++){
            bool bonded = bondPred(vectors[ni1], vectors[ni2]);
            _coordinationStructures[coordType].neighborArray.setNeighborBond(ni1, ni2, bonded);
        }
        _coordinationStructures[coordType].cnaSignatures[ni1] = sigFunc(ni1);
    }
    
    _coordinationStructures[coordType].latticeVectors.assign(vectors, vectors + numNeighbors);
}

void CoordinationStructures::initializeLatticeStructure(
    int latticeType, 
    const Vector3* vectors, 
    int totalVectors,
    CoordinationStructure* coordStruct
){
    _latticeStructures[latticeType].latticeVectors.assign(vectors, vectors + totalVectors);
    _latticeStructures[latticeType].coordStructure = coordStruct;
    _latticeStructures[latticeType].maxNeighbors = coordStruct->numNeighbors;
}

void CoordinationStructures::initializeDiamondStructure(int coordType, int latticeType, const Vector3* vectors, int numNeighbors, int totalVectors){
    _coordinationStructures[coordType].numNeighbors = numNeighbors;
    for(int ni1 = 0; ni1 < numNeighbors; ++ni1){
        _coordinationStructures[coordType].neighborArray.setNeighborBond(ni1, ni1, false);
        double cutoff = (ni1 < 4) ? (sqrt(3.0)*0.25+sqrt(0.5))/2.0 : (1.0+sqrt(0.5))/2.0;

        for(int ni2 = 0; ni2 < 4; ++ni2){
            if(ni1 < 4 && ni2 < 4) _coordinationStructures[coordType].neighborArray.setNeighborBond(ni1, ni2, false);
        }

        for(int ni2 = std::max(ni1 + 1, 4); ni2 < numNeighbors; ++ni2){
            bool bonded = (vectors[ni1] - vectors[ni2]).length() < cutoff;
            _coordinationStructures[coordType].neighborArray.setNeighborBond(ni1, ni2, bonded);
        }

        if(coordType == COORD_HEX_DIAMOND){
            _coordinationStructures[coordType].cnaSignatures[ni1] = (ni1 < 4) ? 0 : ((vectors[ni1].z() == 0) ? 2 : 1);
        } else {
            _coordinationStructures[coordType].cnaSignatures[ni1] = (ni1 < 4) ? 0 : 1;
        }
    }
    _coordinationStructures[coordType].latticeVectors.assign(vectors, vectors + numNeighbors);
    initializeLatticeStructure(latticeType, vectors, totalVectors, &_coordinationStructures[coordType]);
}

void CoordinationStructures::initializeCommonNeighbors(){
	for(auto coordStruct = std::begin(_coordinationStructures); coordStruct != std::end(_coordinationStructures); ++coordStruct){
		for(int neighborIndex = 0; neighborIndex < coordStruct->numNeighbors; neighborIndex++){
			findCommonNeighborsForBond(*coordStruct, neighborIndex);
		}
	}
}

void CoordinationStructures::findCommonNeighborsForBond(CoordinationStructure& coordStruct, int neighborIndex){
	Matrix3 tm;
	tm.column(0) = coordStruct.latticeVectors[neighborIndex];
    bool found = false;
    
    coordStruct.commonNeighbors[neighborIndex][0] = -1;
    coordStruct.commonNeighbors[neighborIndex][1] = -1;
    
    // Special case for SC (Simple Cubic) structure:
    // SC has 6 neighbors along x, y, z. These neighbors don't share common bonds
    // in the traditional CNA sense (each neighbor is perpendicular to 4 others and 
    // opposite to 1). We need any 2 orthogonal vectors that form a complete basis 
    // with the target neighbor vector.
    if(coordStruct.numNeighbors == 6){
        // SC_VECTORS ordering: {+x, -x, +y, -y, +z, -z} -> indices {0,1,2,3,4,5}
        // XOR with 1 gives the opposite: 0↔1, 2↔3, 4↔5
        for(int i1 = 0; i1 < 6 && !found; i1++){
            // Skip the target neighbor and its opposite
            if(i1 == neighborIndex || i1 == (neighborIndex ^ 1)) continue;
            tm.column(1) = coordStruct.latticeVectors[i1];
            
            for(int i2 = i1 + 1; i2 < 6; i2++){
                // Skip target, its opposite, and opposite of i1
                if(i2 == neighborIndex || i2 == (neighborIndex ^ 1)) continue;
                if(i2 == (i1 ^ 1)) continue;
                tm.column(2) = coordStruct.latticeVectors[i2];
                
                if(std::abs(tm.determinant()) > EPSILON){
                    coordStruct.commonNeighbors[neighborIndex][0] = i1;
                    coordStruct.commonNeighbors[neighborIndex][1] = i2;
                    found = true;
                    break;
                }
            }
        }
        if(found) return;
    }
    
    for(int i1 = 0; i1 < coordStruct.numNeighbors && !found; i1++){
        if(!coordStruct.neighborArray.neighborBond(neighborIndex, i1)) continue;
        tm.column(1) = coordStruct.latticeVectors[i1];
        
        for(int i2 = i1 + 1; i2 < coordStruct.numNeighbors; i2++){
            if(!coordStruct.neighborArray.neighborBond(neighborIndex, i2)) continue;
            tm.column(2) = coordStruct.latticeVectors[i2];
            
            if(std::abs(tm.determinant()) > EPSILON){
                coordStruct.commonNeighbors[neighborIndex][0] = i1;
                coordStruct.commonNeighbors[neighborIndex][1] = i2;
                found = true;
                break;
            }
        }
    }
    //assert(found);
}

void CoordinationStructures::initializeSymmetryInformation(){
    for(auto latticeStruct = std::begin(_latticeStructures); 
        latticeStruct != std::end(_latticeStructures); ++latticeStruct){
        
        if(latticeStruct->latticeVectors.empty()) continue;
        
        latticeStruct->primitiveCellInverse = latticeStruct->primitiveCell.inverse();
        generateSymmetryPermutations(*latticeStruct);
        calculateSymmetryProducts(*latticeStruct);
    }
}

void CoordinationStructures::generateSymmetryPermutations(LatticeStructure& latticeStruct){
    const CoordinationStructure& coordStruct = *latticeStruct.coordStructure;

    AnalysisSymmetryUtils::generateSymmetryPermutations(
        coordStruct.latticeVectors,
        coordStruct.numNeighbors,
        latticeStruct.latticeVectors,
        latticeStruct.permutations
    );
}

void CoordinationStructures::findNonCoplanarVectors(const CoordinationStructure& coordStruct, int nindices[3], Matrix3& tm1){
    AnalysisSymmetryUtils::findNonCoplanarVectors(coordStruct.latticeVectors, coordStruct.numNeighbors, nindices, tm1);
}

void CoordinationStructures::calculateSymmetryProducts(LatticeStructure& latticeStruct){
    AnalysisSymmetryUtils::calculateSymmetryProducts(latticeStruct.permutations);
}

void CoordinationStructures::initializeStructures(){
	initializeOther();
	initializeFCC();
	initializeHCP();
	initializeBCC();
    initializeSC();
	initializeCubicDiamond();
	initializeHexagonalDiamond();
	initializeCommonNeighbors();
	initializeSymmetryInformation();
}

}
