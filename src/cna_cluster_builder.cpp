#include <volt/cna_cluster_builder.h>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

namespace Volt{

CNAClusterBuilder::CNAClusterBuilder(
    StructureAnalysis& sa,
    AnalysisContext& context
): ClusterBuilder(sa, context){}

void CNAClusterBuilder::grow(
    Cluster* cluster,
    std::deque<int>& atomsToVisit,
    Matrix_3<double>& orientationV,
    Matrix_3<double>& orientationW,
    int structureType
){
    const int coordinationNumber = _sa.coordinationNumber(structureType);

    while(!atomsToVisit.empty()){
        int currentAtomIndex = atomsToVisit.front();
        atomsToVisit.pop_front();

        int symmetryPermutationIndex = _context.atomSymmetryPermutations->getInt(static_cast<size_t>(currentAtomIndex));

        for(int neighborIndex = 0; neighborIndex < coordinationNumber; neighborIndex++){
            int neighborAtomIndex = _sa.getNeighbor(currentAtomIndex, neighborIndex);

            const Vector3& latticeVector = _sa.latticeVector(
                structureType,
                _sa.symmetryPermutationEntry(structureType, symmetryPermutationIndex, neighborIndex)
            );
            const Vector3& spatialVector = _context.simCell.wrapVector(
                _context.positions->getPoint3(static_cast<size_t>(neighborAtomIndex)) - _context.positions->getPoint3(static_cast<size_t>(currentAtomIndex))
            );

            for(size_t i = 0; i < 3; i++){
                for(size_t j = 0; j < 3; j++){
                    orientationV(i, j) += (latticeVector[j] * latticeVector[i]);
                    orientationW(i, j) += (latticeVector[j] * spatialVector[i]);
                }
            }

            if(_context.atomClusters->getInt(neighborAtomIndex) != 0) continue;
            if(_context.structureTypes->getInt(neighborAtomIndex) != structureType) continue;

            Matrix3 tm1, tm2;
            bool properOverlap = true;

            for(int i = 0; i < 3; i++){
                int atomIndex;
                if(i != 2){
                    int commonNeighbor = _sa.commonNeighborIndex(structureType, neighborIndex, i);
                    atomIndex = _sa.getNeighbor(currentAtomIndex, commonNeighbor);
                    tm1.column(i) = _sa.latticeVector(
                        structureType,
                        _sa.symmetryPermutationEntry(structureType, symmetryPermutationIndex, commonNeighbor)
                    ) - _sa.latticeVector(
                        structureType,
                        _sa.symmetryPermutationEntry(structureType, symmetryPermutationIndex, neighborIndex)
                    );
                }else{
                    atomIndex = currentAtomIndex;
                    tm1.column(i) = -_sa.latticeVector(
                        structureType,
                        _sa.symmetryPermutationEntry(structureType, symmetryPermutationIndex, neighborIndex)
                    );
                }

                int j = _sa.findNeighbor(neighborAtomIndex, atomIndex);
                if(j == -1){
                    properOverlap = false;
                    break;
                }
                tm2.column(i) = _sa.latticeVector(structureType, j);
            }

            if(!properOverlap) continue;

            Matrix3 tm2inverse;
            if(!tm2.inverse(tm2inverse)) continue;

            Matrix3 transition = tm1 * tm2inverse;

            for(int i = 0; i < _sa.symmetryPermutationCount(structureType); i++){
                if(transition.equals(_sa.symmetryTransformation(structureType, i), CA_TRANSITION_MATRIX_EPSILON)){
                    _context.atomClusters->setInt(neighborAtomIndex, cluster->id);
                    cluster->atomCount++;
                    _context.atomSymmetryPermutations->setInt(neighborAtomIndex, i);
                    atomsToVisit.push_back(neighborAtomIndex);
                    break;
                }
            }
        }
    }
}

void CNAClusterBuilder::applyPreferredOrientation(Cluster* cluster){
    double smallestDeviation = std::numeric_limits<double>::max();
    Matrix3 oldOrientation = cluster->orientation;

    for(int symIndex = 0; symIndex < _sa.symmetryPermutationCount(cluster->structure); ++symIndex){
        Matrix3 newOrientation = oldOrientation * _sa.symmetryTransformation(cluster->structure, symIndex).inverse();
        double scaling = std::pow(std::abs(newOrientation.determinant()), 1.0 / 3.0);

        for(const auto& preferredOrientation : _context.preferredCrystalOrientations){
            double deviation = 0;
            for(size_t i = 0; i < 3; i++){
                for(size_t j = 0; j < 3; j++){
                    deviation += std::abs(newOrientation(i, j) / scaling - preferredOrientation(i, j));
                }
            }
            if(deviation < smallestDeviation){
                smallestDeviation = deviation;
                cluster->symmetryTransformation = symIndex;
                cluster->orientation = newOrientation;
            }
        }
    }
}

void CNAClusterBuilder::build(bool dissolveSmallClusters){
    for(size_t seedAtomIndex = 0; seedAtomIndex < _context.atomCount(); seedAtomIndex++){
        if(alreadyProcessedAtom(seedAtomIndex)) continue;

        int structureType = _context.structureTypes->getInt(seedAtomIndex);
        Cluster* cluster = startNew(seedAtomIndex, structureType);

        Matrix_3<double> orientationV = Matrix_3<double>::Zero();
        Matrix_3<double> orientationW = Matrix_3<double>::Zero();
        std::deque<int> atomsToVisit(1, seedAtomIndex);

        grow(cluster, atomsToVisit, orientationV, orientationW, structureType);
        cluster->orientation = Matrix3(orientationW * orientationV.inverse());

        if(structureType == _context.inputCrystalType && !_context.preferredCrystalOrientations.empty()){
            applyPreferredOrientation(cluster);
        }
    }
    
    reorientAtomsToAlign();
    connectClusters();
    formSuperClusters();
    if(dissolveSmallClusters){
        ClusterBuilder::dissolveSmallClusters();
    }
}

bool CNAClusterBuilder::calculateMisorientation(
    int atomIndex,
    int neighbor,
    int neighborIndex,
    Matrix3& outTransition
){
    int structureType = _context.structureTypes->getInt(atomIndex);
    int symIndex = _context.atomSymmetryPermutations->getInt(atomIndex);

    Matrix3 tm1, tm2;
    for(int i = 0; i < 3; i++){
        int ai;
        if(i != 2){
            int cnIdx = _sa.commonNeighborIndex(structureType, neighborIndex, i);
            if(cnIdx < 0){
                return false;
            }

            ai = _sa.getNeighbor(atomIndex, cnIdx);
            tm1.column(i) = _sa.latticeVector(structureType, _sa.symmetryPermutationEntry(structureType, symIndex, cnIdx)) -
                            _sa.latticeVector(structureType, _sa.symmetryPermutationEntry(structureType, symIndex, neighborIndex));
        }else{
            ai = atomIndex;
            tm1.column(i) = -_sa.latticeVector(structureType, _sa.symmetryPermutationEntry(structureType, symIndex, neighborIndex));
        }

        if(_sa.numberOfNeighbors(neighbor) != _sa.coordinationNumber(structureType)){
            return false;
        }

        int j = _sa.findNeighbor(neighbor, ai);
        if(j == -1){
            return false;
        }

        int neighborStructureType = _context.structureTypes->getInt(neighbor);
        int neighborSymIndex = _context.atomSymmetryPermutations->getInt(neighbor);
        tm2.column(i) = _sa.latticeVector(
            neighborStructureType,
            _sa.symmetryPermutationEntry(neighborStructureType, neighborSymIndex, j)
        );
    }

    if(std::abs(tm1.determinant()) < EPSILON){
        return false;
    }

    Matrix3 tm1inv;
    if(!tm1.inverse(tm1inv)){
        return false;
    }

    outTransition = tm2 * tm1inv;
    return true;
}

void CNAClusterBuilder::reorientAtomsToAlign(){
    tbb::parallel_for(tbb::blocked_range<size_t>(0, _context.atomCount()), [this](const tbb::blocked_range<size_t>& r){
        for(size_t atomIndex = r.begin(); atomIndex != r.end(); ++atomIndex){
            int clusterId = _context.atomClusters->getInt(atomIndex);
            if(clusterId == 0){
                continue;
            }

            Cluster* cluster = _sa.clusterGraph().findCluster(clusterId);
            assert(cluster);
            if(cluster->symmetryTransformation == 0){
                continue;
            }

            int oldSymmetry = _context.atomSymmetryPermutations->getInt(atomIndex);
            int newSymmetry = _sa.symmetryInverseProduct(
                cluster->structure,
                oldSymmetry,
                cluster->symmetryTransformation
            );

            _context.atomSymmetryPermutations->setInt(atomIndex, newSymmetry);
        }
    });
}

void CNAClusterBuilder::connectClusters(){
    std::vector<std::vector<int>> extras(_context.atomCount());

    for(size_t atomIndex = 0; atomIndex < _context.atomCount(); ++atomIndex){
        int clusterId = _context.atomClusters->getInt(atomIndex);
        if(clusterId == 0){
            continue;
        }

        Cluster* cluster1 = _sa.clusterGraph().findCluster(clusterId);
        const int nn = _sa.numberOfNeighbors(atomIndex);

        for(int ni = 0; ni < nn; ++ni){
            int neighbor = _sa.getNeighbor(atomIndex, ni);
            if(neighbor < 0 || neighbor == static_cast<int>(atomIndex)){
                continue;
            }

            int neighborClusterId = _context.atomClusters->getInt(neighbor);
            if(neighborClusterId == 0){
                extras[neighbor].push_back(static_cast<int>(atomIndex));
                continue;
            }

            if(neighborClusterId == cluster1->id){
                continue;
            }

            Cluster* cluster2 = _sa.clusterGraph().findCluster(neighborClusterId);
            if(ClusterTransition* existing = cluster1->findTransition(cluster2)){
                existing->area++;
                existing->reverse->area++;
                continue;
            }

            Matrix3 transition;
            if(!calculateMisorientation(static_cast<int>(atomIndex), neighbor, ni, transition)){
                continue;
            }
            if(!transition.isOrthogonalMatrix()){
                continue;
            }

            ClusterTransition* transitionLink = _sa.clusterGraph().createClusterTransition(cluster1, cluster2, transition);
            transitionLink->area++;
            transitionLink->reverse->area++;
        }
    }

    _sa.appendNeighbors(extras);
    spdlog::info("Number of cluster transitions: {}", _sa.clusterGraph().clusterTransitions().size());
}

void CNAClusterBuilder::processDefectCluster(Cluster* defectCluster){
    for(ClusterTransition* transition = defectCluster->transitions; transition; transition = transition->next){
        if(transition->cluster2->structure != _context.inputCrystalType || transition->distance != 1){
            continue;
        }
        for(ClusterTransition* sibling = transition->next; sibling; sibling = sibling->next) {
            if(sibling->cluster2->structure != _context.inputCrystalType || sibling->distance != 1){
                continue;
            }
            if(sibling->cluster2 == transition->cluster2){
                continue;
            }

            Matrix3 misorientation = sibling->tm * transition->reverse->tm;

            for(int symIndex = 0; symIndex < _sa.symmetryPermutationCount(sibling->cluster2->structure); ++symIndex){
                if(_sa.symmetryTransformation(sibling->cluster2->structure, symIndex).equals(misorientation, 1e-6)){
                    _sa.clusterGraph().createClusterTransition(transition->cluster2, sibling->cluster2, misorientation, 2);
                    break;
                }
            }
        }
    }
}

void CNAClusterBuilder::formSuperClusters(){
    const size_t oldTransitionCount = _sa.clusterGraph().clusterTransitions().size();

    for(Cluster* cluster : _sa.clusterGraph().clusters()){
        if(!cluster || cluster->id == 0){
            continue;
        }
        cluster->rank = 0;
    }

    for(Cluster* cluster : _sa.clusterGraph().clusters()){
        if(!cluster || cluster->id == 0){
            continue;
        }
        if(cluster->structure != _context.inputCrystalType){
            processDefectCluster(cluster);
        }
    }

    const size_t newTransitionCount = _sa.clusterGraph().clusterTransitions().size();

    for(size_t i = oldTransitionCount; i < newTransitionCount; i++){
        ClusterTransition* transition = _sa.clusterGraph().clusterTransitions()[i];

        Cluster* parent1 = getParentGrain(transition->cluster1);
        Cluster* parent2 = getParentGrain(transition->cluster2);

        if(parent1 == parent2){
            continue;
        }

        ClusterTransition* parentTransition = transition;

        if(parent2 != transition->cluster2){
            parentTransition = _sa.clusterGraph().concatenateClusterTransitions(
                parentTransition,
                transition->cluster2->parentTransition
            );
        }

        if(parent1 != transition->cluster1){
            parentTransition = _sa.clusterGraph().concatenateClusterTransitions(
                transition->cluster1->parentTransition->reverse,
                parentTransition
            );
        }

        if(parent1->rank > parent2->rank){
            parent2->parentTransition = parentTransition->reverse;
            continue;
        }

        parent1->parentTransition = parentTransition;

        if(parent1->rank == parent2->rank){
            parent2->rank++;
        }
    }

    for(Cluster* cluster : _sa.clusterGraph().clusters()){
        getParentGrain(cluster);
    }
}

}
