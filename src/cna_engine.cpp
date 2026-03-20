#include <volt/cna_engine.h>
#include <volt/core/neighbor_ordering.h>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <cstdint>

namespace Volt{

CommonNeighborAnalysisEngine::CommonNeighborAnalysisEngine(AnalysisContext& context, bool identifyPlanarDefects)
    : _context(context)
    , _identifyPlanarDefects(identifyPlanarDefects){}

int CommonNeighborAnalysisEngine::coordinationNumber() const{
    switch(_context.inputCrystalType){
        case LATTICE_FCC:
        case LATTICE_HCP:
            return 12;
        case LATTICE_BCC:
            return 14;
        case LATTICE_CUBIC_DIAMOND:
        case LATTICE_HEX_DIAMOND:
            return 16;
        case LATTICE_SC:
            return 6;
        default:
            return 0;
    }
}

double CommonNeighborAnalysisEngine::computeLocalCutoff(
    const NearestNeighborFinder& neighList,
    const NearestNeighborFinder::Query<MAX_NEIGHBORS>& neighQuery,
    int numNeighbors,
    int coordinationCount,
    int particleIndex,
    int* neighborIndices,
    Vector3* neighborVectors,
    NeighborBondArray& neighborArray
) const{
    double localScaling = 0.0;
    double localCutoff = 0.0;

    switch(_context.inputCrystalType){
        case LATTICE_FCC:
        case LATTICE_HCP:
            for(int neighbor = 0; neighbor < 12; neighbor++){
                localScaling += sqrt(neighQuery.results()[neighbor].distanceSq);
            }
            localScaling /= 12;
            localCutoff = localScaling * (1.0f + sqrt(2.0f)) * 0.5f;
            break;
        case LATTICE_BCC:
            for(int neighbor = 0; neighbor < 8; neighbor++){
                localScaling += sqrt(neighQuery.results()[neighbor].distanceSq);
            }
            localScaling /= 8;
            localCutoff = localScaling / (sqrt(3.0) / 2.0) * 0.5 * (1.0 + sqrt(2.0));
            break;
        case LATTICE_CUBIC_DIAMOND:
        case LATTICE_HEX_DIAMOND: {
            int outputIndex = 4;
            for(int i = 0; i < 4; i++){
                const Vector3& v0 = neighQuery.results()[i].delta;
                neighborVectors[i] = v0;
                neighborIndices[i] = neighQuery.results()[i].index;

                NearestNeighborFinder::Query<MAX_NEIGHBORS> neighQuery2(neighList);
                neighQuery2.findNeighbors(neighList.particlePos(neighborIndices[i]));
                if(neighQuery2.results().size() < 4) return 0.0;

                for(int j = 0; j < 4; j++){
                    Vector3 v = v0 + neighQuery2.results()[j].delta;
                    if(neighQuery2.results()[j].index == particleIndex && v.isZero()) continue;
                    if(outputIndex == 16) return 0.0;

                    neighborIndices[outputIndex] = neighQuery2.results()[j].index;
                    neighborVectors[outputIndex] = v;
                    neighborArray.setNeighborBond(i, outputIndex, true);
                    outputIndex++;
                }

                if(outputIndex != (i * 3) + 7) return 0.0;
            }

            localScaling = 0.0;
            for(int neighbor = 4; neighbor < 16; neighbor++){
                localScaling += neighborVectors[neighbor].length();
            }

            localScaling /= 12;
            localCutoff = localScaling * 1.2071068;
            break;
        }
        default:
            return 0.0;
    }

    const double localCutoffSquared = localCutoff * localCutoff;

    switch(_context.inputCrystalType){
        case LATTICE_FCC:
        case LATTICE_HCP:
        case LATTICE_BCC:
        case LATTICE_SC:
            if(numNeighbors > coordinationCount && neighQuery.results()[coordinationCount].distanceSq <= localCutoffSquared){
                return 0.0;
            }

            for(int ni1 = 0; ni1 < coordinationCount; ni1++){
                neighborIndices[ni1] = neighQuery.results()[ni1].index;
                neighborVectors[ni1] = neighQuery.results()[ni1].delta;
                neighborArray.setNeighborBond(ni1, ni1, false);
                for(int ni2 = ni1 + 1; ni2 < coordinationCount; ni2++){
                    neighborArray.setNeighborBond(
                        ni1,
                        ni2,
                        (neighQuery.results()[ni1].delta - neighQuery.results()[ni2].delta).squaredLength() <= localCutoffSquared
                    );
                }
            }
            break;
        case LATTICE_CUBIC_DIAMOND:
        case LATTICE_HEX_DIAMOND:
            for(int ni1 = 4; ni1 < coordinationCount; ni1++){
                for(int ni2 = ni1 + 1; ni2 < coordinationCount; ni2++){
                    auto distance = (neighborVectors[ni1] - neighborVectors[ni2]);
                    bool isBonded = distance.squaredLength() <= localCutoffSquared;
                    neighborArray.setNeighborBond(ni1, ni2, isBonded);
                }
            }
            break;
        default:
            return 0.0;
    }

    return localCutoff;
}

double CommonNeighborAnalysisEngine::determineLocalStructure(
    const NearestNeighborFinder& neighList,
    int particleIndex
) const{
    std::array<int, MAX_NEIGHBORS> neighborIndices;
    std::array<Vector3, MAX_NEIGHBORS> neighborVectors;
    std::array<int, MAX_NEIGHBORS> cnaSignatures;
    std::array<int, MAX_NEIGHBORS> neighborMapping;
    std::array<int, MAX_NEIGHBORS> previousMapping;
    NeighborBondArray neighborArray;

    NearestNeighborFinder::Query<MAX_NEIGHBORS> neighQuery(neighList);
    neighQuery.findNeighbors(neighList.particlePos(particleIndex));
    const int numNeighbors = static_cast<int>(neighQuery.results().size());
    const int coordinationCount = coordinationNumber();
    if(numNeighbors < coordinationCount){
        return 0.0;
    }

    const double localCutoff = computeLocalCutoff(
        neighList,
        neighQuery,
        numNeighbors,
        coordinationCount,
        particleIndex,
        neighborIndices.data(),
        neighborVectors.data(),
        neighborArray
    );
    if(localCutoff == 0.0){
        return 0.0;
    }

    for(int n = 0; n < coordinationCount; ++n){
        neighborMapping[n] = n;
        previousMapping[n] = -1;
    }

    CoordinationStructureType coordinationType = CommonNeighborAnalysis::computeCoordinationType(
        neighborArray,
        coordinationCount,
        cnaSignatures.data(),
        _context.inputCrystalType,
        _identifyPlanarDefects
    );
    if(coordinationType == COORD_OTHER){
        return 0.0;
    }

    bool found = CommonNeighborAnalysis::findMatchingNeighborPermutation(
        coordinationType,
        neighborMapping.data(),
        previousMapping.data(),
        coordinationCount,
        cnaSignatures.data(),
        neighborArray,
        CoordinationStructures::_coordinationStructures
    );
    if(!found){
        return 0.0;
    }

    StructureType atomStructure = StructureType::OTHER;
    switch(coordinationType){
        case COORD_CUBIC_DIAMOND:
            atomStructure = StructureType::CUBIC_DIAMOND;
            break;
        case COORD_HEX_DIAMOND:
            atomStructure = StructureType::HEX_DIAMOND;
            break;
        case COORD_FCC:
            atomStructure = StructureType::FCC;
            break;
        case COORD_HCP:
            atomStructure = StructureType::HCP;
            break;
        case COORD_BCC:
            atomStructure = StructureType::BCC;
            break;
        case COORD_SC:
            atomStructure = StructureType::SC;
            break;
        default:
            atomStructure = StructureType::OTHER;
            break;
    }

    _context.structureTypes->setInt(particleIndex, static_cast<int>(atomStructure));
    storeNeighborOrdering(
        neighborMapping,
        coordinationCount,
        atomStructure,
        particleIndex
    );

    for(int i = 0; i < coordinationCount; i++){
        const Vector3& neighborVector = neighborVectors[neighborMapping[i]];
        for(int dim = 0; dim < 3; dim++){
            if(_context.simCell.pbcFlags()[dim]){
                if(std::abs(_context.simCell.inverseMatrix().prodrow(neighborVector, dim)) >= 0.5 + EPSILON){
                    CoordinationStructures::generateCellTooSmallError(dim);
                }
            }
        }
    }

    return localCutoff;
}

void CommonNeighborAnalysisEngine::storeNeighborOrdering(
    const std::array<int, MAX_NEIGHBORS>& neighborMapping,
    int coordinationCount,
    StructureType atomStructure,
    int particleIndex
) const{
    if(!_context.correspondences || coordinationCount <= 0){
        return;
    }

    if(!describeNeighborOrdering(atomStructure).isSingleShell()){
        return;
    }

    NeighborOrderingStorage encodedOrdering{};
    encodedOrdering[0] = 0;
    for(int i = 0; i < coordinationCount; ++i){
        encodedOrdering[static_cast<size_t>(i + 1)] = static_cast<int8_t>(neighborMapping[static_cast<size_t>(i)] + 1);
    }

    std::uint64_t orderingCode = 0;
    if(!encodeNeighborOrdering(atomStructure, coordinationCount, encodedOrdering, orderingCode)){
        return;
    }

    _context.correspondences->setInt64(static_cast<size_t>(particleIndex), static_cast<std::int64_t>(orderingCode));
}

void CommonNeighborAnalysisEngine::identifyStructures(){
    static std::once_flag init_flag;
    std::call_once(init_flag, [](){
        CoordinationStructures::initializeStructures();
    });

    std::fill(
        _context.structureTypes->dataInt(),
        _context.structureTypes->dataInt() + _context.structureTypes->size(),
        static_cast<int>(StructureType::OTHER)
    );

    NearestNeighborFinder neighFinder(MAX_NEIGHBORS);
    if(!neighFinder.prepare(_context.positions, _context.simCell, _context.particleSelection)){
        throw std::runtime_error("Error in neighFinder.prepare(...)");
    }

    const size_t N = _context.atomCount();
    tbb::parallel_for(tbb::blocked_range<size_t>(0, N), [&](const tbb::blocked_range<size_t>& r){
        for(size_t index = r.begin(); index != r.end(); ++index){
            determineLocalStructure(neighFinder, static_cast<int>(index));
        }
    });

    invalidateStatistics();
}

void CommonNeighborAnalysisEngine::calculateStructureStatistics() const{
    _structureStatistics.clear();
    const size_t N = _context.atomCount();
    for(size_t i = 0; i < N; ++i){
        int structureType = _context.structureTypes->getInt(i);
        _structureStatistics[structureType]++;
    }
    _statisticsValid = true;
}

void CommonNeighborAnalysisEngine::invalidateStatistics(){
    _statisticsValid = false;
}

std::string CommonNeighborAnalysisEngine::getStructureTypeName(int structureType) const{
    return std::string(structureTypeName(structureType));
}

json CommonNeighborAnalysisEngine::buildMainListing() const{
    if(!_statisticsValid){
        calculateStructureStatistics();
    }

    const int N = static_cast<int>(_context.atomCount());
    const double invN = (N > 0) ? (100.0 / static_cast<double>(N)) : 0.0;

    int totalIdentified = 0;
    int unidentified = 0;
    auto itOther = _structureStatistics.find(static_cast<int>(StructureType::OTHER));
    if(itOther != _structureStatistics.end()){
        unidentified = itOther->second;
    }

    json mainListing = json::object();
    mainListing["total_atoms"] = N;
    mainListing["analysis_method"] = "CNA";

    for(const auto& [structureType, count] : _structureStatistics){
        std::string name = getStructureTypeName(structureType);
        std::transform(name.begin(), name.end(), name.begin(), ::tolower);
        mainListing[name + "_count"] = count;
        mainListing[name + "_percentage"] = static_cast<double>(count) * invN;
        if(structureType != static_cast<int>(StructureType::OTHER) &&
           structureType != static_cast<int>(CoordinationStructureType::COORD_OTHER)){
            totalIdentified += count;
        }
    }

    mainListing["total_identified"] = totalIdentified;
    mainListing["total_unidentified"] = unidentified;
    mainListing["identification_rate"] = static_cast<double>(totalIdentified) * invN;
    mainListing["unique_structure_types"] = static_cast<int>(_structureStatistics.size());

    return mainListing;
}

json CommonNeighborAnalysisEngine::getPerAtomProperties(const LammpsParser::Frame& frame) const{
    json perAtom = json::array();

    for(size_t i = 0; i < static_cast<size_t>(frame.natoms); ++i){
        const int structureType = _context.structureTypes->getInt(i);

        json atom;
        atom["id"] = i < frame.ids.size() ? frame.ids[i] : static_cast<int>(i);
        atom["structure_type"] = structureType;
        atom["structure_name"] = getStructureTypeName(structureType);

        if(i < frame.positions.size()){
            const auto& pos = frame.positions[i];
            atom["pos"] = {pos.x(), pos.y(), pos.z()};
        }else{
            atom["pos"] = {0.0, 0.0, 0.0};
        }

        perAtom.push_back(atom);
    }

    return perAtom;
}

}
