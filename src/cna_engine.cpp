#include <volt/cna_engine.h>
#include <volt/core/neighbor_ordering.h>
#include <volt/coordination_structures_utils.h>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <cstdint>

namespace Volt{

CommonNeighborAnalysisEngine::CommonNeighborAnalysisEngine(AnalysisContext& context, bool identifyPlanarDefects)
    : _context(context), _identifyPlanarDefects(identifyPlanarDefects){}

void CommonNeighborAnalysisEngine::perform(){
    identifyStructures();
    calculateStructureStatistics();
}

double CommonNeighborAnalysisEngine::determineLocalStructure(
    const NearestNeighborFinder& neighList,
    int particleIndex
) const{
    CnaLocalStructureUtils::LocalStructureMatch match;
    if(!CnaLocalStructureUtils::determineLocalStructure(
        _context.inputCrystalType,
        _identifyPlanarDefects,
        neighList,
        particleIndex,
        CoordinationStructures::_coordinationStructures,
        match
    )){
        return 0.0;
    }

    _context.structureTypes->setInt(particleIndex, static_cast<int>(match.structureType));
    storeNeighborOrdering(
        match.neighborMapping,
        match.coordinationNumber,
        match.structureType,
        particleIndex
    );

    return match.localCutoff;
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
    ensureCoordinationStructuresInitialized();

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
}

void CommonNeighborAnalysisEngine::calculateStructureStatistics() const{
}

std::string CommonNeighborAnalysisEngine::getStructureTypeName(int structureType) const{
    return std::string(structureTypeName(structureType));
}

json CommonNeighborAnalysisEngine::buildMainListing() const{
    return json::object();
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
