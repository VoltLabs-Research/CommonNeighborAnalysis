#include <volt/analysis/cna_classifier.h>
#include <volt/structures/crystal_topology_registry.h>

#include <algorithm>
#include <cassert>
#include <cstring>
#include <vector>

namespace Volt{

struct CnaSignatureDescriptor{
    int commonNeighborCount = 0;
    int neighborBondCount = 0;
    int maxChainLength = 0;
};

bool operator==(const CnaSignatureDescriptor& lhs, const CnaSignatureDescriptor& rhs){
    return lhs.commonNeighborCount == rhs.commonNeighborCount &&
        lhs.neighborBondCount == rhs.neighborBondCount &&
        lhs.maxChainLength == rhs.maxChainLength;
}

CnaSignatureDescriptor buildSignatureDescriptor(
    const NeighborBondArray& neighborArray,
    int neighborIndex,
    int coordinationNumber
){
    CnaSignatureDescriptor descriptor;
    unsigned int commonNeighbors = 0;
    descriptor.commonNeighborCount = CommonNeighborAnalysis::findCommonNeighbors(
        neighborArray,
        neighborIndex,
        commonNeighbors,
        coordinationNumber
    );

    CNAPairBond neighborBonds[MAX_NEIGHBORS * MAX_NEIGHBORS];
    descriptor.neighborBondCount = CommonNeighborAnalysis::findNeighborBonds(
        neighborArray,
        commonNeighbors,
        coordinationNumber,
        neighborBonds
    );
    descriptor.maxChainLength = CommonNeighborAnalysis::calcMaxChainLength(
        neighborBonds,
        descriptor.neighborBondCount
    );
    return descriptor;
}

NeighborBondArray makeExpectedNeighborArray(const CrystalTopologyEntry& topology){
    NeighborBondArray neighborArray;
    for(int index = 0; index < topology.coordinationNumber; ++index){
        neighborArray.neighborArray[index] = topology.neighborBondRows[static_cast<std::size_t>(index)];
    }
    return neighborArray;
}

bool tryMapDescriptorCode(
    const std::vector<std::pair<CnaSignatureDescriptor, int>>& descriptorToCode,
    const CnaSignatureDescriptor& descriptor,
    int& code
){
    for(const auto& pair : descriptorToCode){
        if(pair.first == descriptor){
            code = pair.second;
            return true;
        }
    }
    return false;
}

bool tryMatchCandidateTopology(
    const CrystalTopologyEntry& topology,
    const std::vector<CnaSignatureDescriptor>& actualDescriptors,
    int* cnaSignatures
){
    if(static_cast<int>(topology.cnaSignatureCodes.size()) != topology.coordinationNumber){
        return false;
    }

    const NeighborBondArray expectedNeighborArray = makeExpectedNeighborArray(topology);
    std::vector<std::pair<CnaSignatureDescriptor, int>> descriptorToCode;
    descriptorToCode.reserve(static_cast<std::size_t>(topology.coordinationNumber));

    for(int index = 0; index < topology.coordinationNumber; ++index){
        const CnaSignatureDescriptor expectedDescriptor = buildSignatureDescriptor(
            expectedNeighborArray,
            index,
            topology.coordinationNumber
        );
        const int expectedCode = topology.cnaSignatureCodes[static_cast<std::size_t>(index)];

        bool found = false;
        for(const auto& pair : descriptorToCode){
            if(pair.first == expectedDescriptor){
                if(pair.second != expectedCode){
                    return false;
                }
                found = true;
                break;
            }
        }

        if(!found){
            descriptorToCode.emplace_back(expectedDescriptor, expectedCode);
        }
    }

    std::vector<int> actualCodes(actualDescriptors.size(), -1);
    for(std::size_t index = 0; index < actualDescriptors.size(); ++index){
        if(!tryMapDescriptorCode(descriptorToCode, actualDescriptors[index], actualCodes[index])){
            return false;
        }
    }

    std::vector<int> sortedActualCodes = actualCodes;
    std::vector<int> expectedCodes = topology.cnaSignatureCodes;
    std::sort(sortedActualCodes.begin(), sortedActualCodes.end());
    std::sort(expectedCodes.begin(), expectedCodes.end());
    if(sortedActualCodes != expectedCodes){
        return false;
    }

    std::copy(actualCodes.begin(), actualCodes.end(), cnaSignatures);
    return true;
}

std::vector<const CrystalTopologyEntry*> candidateTopologiesForCna(
    LatticeStructureType inputCrystalType,
    bool identifyPlanarDefects
){
    std::vector<const CrystalTopologyEntry*> candidates;
    const auto* inputTopology = crystalTopologyByLatticeType(static_cast<int>(inputCrystalType));
    if(!inputTopology){
        return candidates;
    }

    candidates.push_back(inputTopology);
    if(!identifyPlanarDefects){
        return candidates;
    }

    for(const auto& entry : crystalTopologyRegistry().entries()){
        if(entry.coordinationType <= 0){
            continue;
        }
        if(entry.coordinationNumber != inputTopology->coordinationNumber){
            continue;
        }
        if(entry.coordinationType == inputTopology->coordinationType){
            continue;
        }
        if(static_cast<int>(entry.cnaSignatureCodes.size()) != entry.coordinationNumber){
            continue;
        }
        candidates.push_back(&entry);
    }

    return candidates;
}

int CommonNeighborAnalysis::findCommonNeighbors(
    const NeighborBondArray& neighborArray,
    int neighborIndex,
    unsigned int &commonNeighbors,
    int
){
    commonNeighbors = neighborArray.neighborArray[neighborIndex];
    return __builtin_popcount(commonNeighbors);
}

bool CommonNeighborAnalysis::findMatchingNeighborPermutation(
    CoordinationStructureType coordinationType,
    int* neighborMapping,
    int* previousMapping,
    int coordinationNumber,
    const int* cnaSignatures,
    const NeighborBondArray& neighborArray,
    const CoordinationStructure* coordinationStructures
){
    const CoordinationStructure& coordStructure = coordinationStructures[coordinationType];

    for(;;){
        int ni1 = 0;

        while(neighborMapping[ni1] == previousMapping[ni1]){
            ni1++;
            assert(ni1 < coordinationNumber);
        }

        for(; ni1 < coordinationNumber; ni1++){
            int atomNeighborIndex1 = neighborMapping[ni1];
            previousMapping[ni1] = atomNeighborIndex1;

            if(cnaSignatures[atomNeighborIndex1] != coordStructure.cnaSignatures[ni1]){
                break;
            }

            int ni2;
            for(ni2 = 0; ni2 < ni1; ni2++){
                int atomNeighborIndex2 = neighborMapping[ni2];
                if(neighborArray.neighborBond(atomNeighborIndex1, atomNeighborIndex2) !=
                   coordStructure.neighborArray.neighborBond(ni1, ni2)){
                    break;
                }
            }

            if(ni2 != ni1){
                break;
            }
        }

        if(ni1 == coordinationNumber){
            return true;
        }

        bitmapSort(neighborMapping + ni1 + 1, neighborMapping + coordinationNumber, coordinationNumber);
        if(!std::next_permutation(neighborMapping, neighborMapping + coordinationNumber)){
            assert(false);
            return false;
        }
    }
}

CoordinationStructureType CommonNeighborAnalysis::computeCoordinationType(
    const NeighborBondArray& neighborArray,
    int coordinationNumber,
    int* cnaSignatures,
    LatticeStructureType inputCrystalType,
    bool identifyPlanarDefects
){
    const std::vector<const CrystalTopologyEntry*> candidates = candidateTopologiesForCna(
        inputCrystalType,
        identifyPlanarDefects
    );
    if(candidates.empty()){
        return COORD_OTHER;
    }

    std::vector<CnaSignatureDescriptor> actualDescriptors(static_cast<std::size_t>(coordinationNumber));
    for(int neighborIndex = 0; neighborIndex < coordinationNumber; ++neighborIndex){
        actualDescriptors[static_cast<std::size_t>(neighborIndex)] = buildSignatureDescriptor(
            neighborArray,
            neighborIndex,
            coordinationNumber
        );
    }

    for(const auto* candidate : candidates){
        if(!candidate || candidate->coordinationNumber != coordinationNumber){
            continue;
        }
        if(tryMatchCandidateTopology(*candidate, actualDescriptors, cnaSignatures)){
            return static_cast<CoordinationStructureType>(candidate->coordinationType);
        }
    }

    return COORD_OTHER;
}

int CommonNeighborAnalysis::findNeighborBonds(
    const NeighborBondArray& neighborArray,
    unsigned int commonNeighbors,
    int numNeighbors,
    CNAPairBond* neighborBonds
){
    int numBonds = 0;
    unsigned int nib[32];
    int nibn = 0;
    unsigned int ni1b = 1;

    for(int ni1 = 0; ni1 < numNeighbors; ni1++, ni1b <<= 1){
        if(commonNeighbors & ni1b){
            unsigned int b = commonNeighbors & neighborArray.neighborArray[ni1];
            for(int n = 0; n < nibn; n++){
                if(b & nib[n]){
                    neighborBonds[numBonds++] = ni1b | nib[n];
                }
            }
            nib[nibn++] = ni1b;
        }
    }

    return numBonds;
}

int CommonNeighborAnalysis::getAdjacentBonds(
    unsigned int atom,
    CNAPairBond* bondsToProcess,
    int& numBonds,
    unsigned int& atomsToProcess,
    unsigned int& atomsProcessed
){
    int adjacentBonds = 0;
    for(int b = numBonds - 1; b >= 0; --b){
        if(atom & bondsToProcess[b]){
            ++adjacentBonds;
            atomsToProcess |= bondsToProcess[b] & (~atomsProcessed);
            memmove(&bondsToProcess[b], &bondsToProcess[b+1], sizeof(CNAPairBond) * (numBonds - b - 1));
            --numBonds;
        }
    }
    return adjacentBonds;
}

int CommonNeighborAnalysis::calcMaxChainLength(CNAPairBond* neighborBonds, int numBonds){
    int maxChainLength = 0;

    while(numBonds){
        numBonds--;
        unsigned int atomsToProcess = neighborBonds[numBonds];
        unsigned int atomsProcessed = 0;
        int clusterSize = 1;
        do{
            int nextAtomIndex = __builtin_ctz(atomsToProcess);
            unsigned int nextAtom = 1 << nextAtomIndex;
            atomsProcessed |= nextAtom;
            atomsToProcess &= ~nextAtom;

            clusterSize += getAdjacentBonds(nextAtom, neighborBonds, numBonds, atomsToProcess, atomsProcessed);
        }while(atomsToProcess);
        if(clusterSize > maxChainLength){
            maxChainLength = clusterSize;
        }
    }

    return maxChainLength;
}

}
