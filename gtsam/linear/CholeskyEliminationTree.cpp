/**
* @file    CholeskyEliminationTree.h
* @brief   Elimination tree structure to perform symbolic factorization and Cholesky factorization
* @author  Roger Hsiao (rogerhh)
* @date    Feb. 8, 2023
*/

// GTSAM TODO: 
// 1. Reorder based on clique structure. I.e. after reordering, group cliques together
// 2. Marginalization
// 3. Ability to add factors to variables in the middle
// 4. After symbolic elimination, allocate contiguous workspace for hardware
//


#include <gtsam/inference/Ordering.h>
#include <gtsam/linear/CholeskyEliminationTree.h>
#include <gtsam/linear/CholeskyEliminationTreeNode.h>
#include <gtsam/linear/CholeskyEliminationTreeClique.h>
#include <gtsam/base/LowerTriangularColumnMatrix.h>
#include <gtsam/3rdparty/CCOLAMD/Include/ccolamd.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <algorithm>
#include <fstream>
#include <cmath>

using namespace std;

namespace gtsam {

using ColMajorMatrix = CholeskyEliminationTree::ColMajorMatrix;
using BlockIndexVector = CholeskyEliminationTree::BlockIndexVector;
const size_t CholeskyEliminationTree::LAST_ROW = -1;

CholeskyEliminationTree::CholeskyEliminationTree() 
    : delta_(0, 1), orderingLess(this) { 
}

void CholeskyEliminationTree::addVariables(const Values& newTheta) {
    // cout << "[CholeskyEliminationTree] addVariables()" << endl;
    for(const auto& keyValPair : newTheta) {
        const Key& key = keyValPair.key;
        const Value& val = keyValPair.value;
        const size_t dim = val.dim();
        assert(key == nodes_.size());
        cholesky_.addColumn(key, dim);
        jacobian_.addColumn(key, dim);
        sharedNode newNode = make_shared<Node>(this, key, ordering_version_);
        nodes_.push_back(newNode);
        // We cannot make sharedNode part of the constructor
        // Because at construct time, Clique is not pointed to by shared_ptr
        sharedClique newClique = make_shared<Clique>(this);
        newClique->addNode(newNode);
        ordering_.push_back(ordering_.size());
        orderingToKey_.push_back(key);
        descendants_.push_back(vector<Key>());
        changedDescendants_.push_back(vector<Key>());
        // descendants_.push_back(vector<pair<Key, size_t>>());
        // changedDescendants_.push_back(vector<pair<Key, size_t>>());
        markedStatus_.push_back(NEW);
        backSolveKeys_.push_back(true);
        is_reordered_.push_back(false);

        bool alloc = delta_.preallocateBlock(key, dim, true);
        assert(alloc);

        // // DEBUG
        // descendantStatus.push_back(unordered_map<Key, DescendantStatus>());
    }
    delta_.resolveAllocate();
} 

void CholeskyEliminationTree::markAffectedKeys(const NonlinearFactorGraph& nonlinearFactors,
                                               const FactorIndices& newFactorIndices,
                                               const KeySet& relinKeys, 
                                               const std::optional<FastList<Key>>& extraKeys,
                                               KeySet* affectedKeys,
                                               KeySet* observedKeys) {
    // cout << "[CholeskyEliminationTree] markAffectedKeys()" << endl;
    affectedKeys->clear();
    observedKeys->clear();

    // cout << "RelinKeys: ";
    // for(const Key relinKey : relinKeys) {
    //     cout << relinKey << " ";
    // }
    // cout << endl;

    // RelinKeys should be processed before we add in factors because we only need to
    // relinearize old factors
    for(const Key relinKey : relinKeys) {
        sharedNode relinNode = nodes_[relinKey];
        relinNode->relinearize = true;

        for(const FactorIndex factorIndex : relinNode->factorIndices) {
            factorLinearizeStatus_[factorIndex] = RELINEARIZED;
        }

        auto it = relinNode->factorColStructure.begin();
        auto itSelf = relinNode->factorColStructure.find(relinKey);
        auto itEnd = relinNode->factorColStructure.end();

        // All keys that interact with this key but are lower than this key
        while(it != itSelf) {
            nodes_[*it]->changedFactorColStructure.insert(relinKey);
            affectedKeys->insert(*it);
            it++;
        }
        // All keys that interact with this key but are higher
        while(it != itEnd) {
            relinNode->changedFactorColStructure.insert(*it);
            affectedKeys->insert(*it);
            it++;
        }
    }

    /*
    for(const Key relinKey : relinKeys) {
        sharedNode relinNode = nodes_[relinKey];
        relinNode->relinearize = true;
        for(const FactorIndex factorIndex : relinNode->factorIndices) {
            if(factorLinearizeStatus_[factorIndex] != LINEARIZED) {
                // If linear status of that factor is not linearized
                // Then it's either a new factor or a factor that was already counted
                continue;
            }
            factorLinearizeStatus_[factorIndex] = RELINEARIZED;
            // cout << "set factor " << factorIndex << " to RELINEARIZED" << endl;
            sharedFactor factor = nonlinearFactors[factorIndex];

            // Sort factor keys first
            auto factorKeys = factor->keys();
            // Can't sort factor->keys() itself because it messes up linearization
            sort(factorKeys.begin(), factorKeys.end(), orderingLess);
            
            // Set up changedFactorColStructure and add to affectedKeys
            // Lower keys affect higher keys but not the other way around
            // This needs to be redone after reordering!
            for(auto it1 = factorKeys.begin(); it1 != factorKeys.end(); it1++) {
                Key lowerKey = *it1;
                affectedKeys->insert(lowerKey); // all keys are affected
                sharedNode lowerNode = nodes_[lowerKey];
                for(auto it2 = it1; it2 != factorKeys.end(); it2++) {
                    Key higherKey = *it2;
                    lowerNode->changedFactorColStructure.push_back(higherKey);
                    cout << "node " << lowerKey << " push back " << higherKey << endl;
                }
            }
        }
    }
    */

    // Data structure used to remove duplicate keys and sort
    unordered_map<Key, set<Key>> sortedFactorKeys;

    for(const FactorIndex factorIndex : newFactorIndices) {
        assert(factorIndex == factorLinearizeStatus_.size());
        factorLinearizeStatus_.push_back(UNLINEARIZED);
        sharedFactor factor = nonlinearFactors[factorIndex];
        factors_.push_back(factor);
        jacobian_.preallocateBlock(factorIndex, -1, factor->dim(), true);

        for(Key k1 : factor->keys()) {
            sharedNode node1 = nodes_[k1];
            node1->factorIndices.push_back(factorIndex);
            affectedKeys->insert(k1);
            observedKeys->insert(k1);
            jacobian_.preallocateBlock(factorIndex, k1, factor->dim(), true);
            for(Key k2 : factor->keys()) {
                node1->factorColStructure.insert(k2);
            }
        }

        // // Sort factor keys first
        // auto factorKeys = factor->keys();
        // // Can't sort factor->keys() itself because it messes up linearization
        // sort(factorKeys.begin(), factorKeys.end(), orderingLess);
        // for(auto it1 = factorKeys.begin(); it1 != factorKeys.end(); it1++) {
        //     Key lowerKey = *it1;
        //     sharedNode lowerNode = nodes_[lowerKey];
        //     lowerNode->factorIndices.push_back(factorIndex);
        //     affectedKeys->insert(lowerKey);
        // cout << "jacobian preallocate " << factorIndex << " " << lowerKey << endl;
        //     jacobian_.preallocateBlock(factorIndex, lowerKey, factor->dim(), true);

        //     for(auto it2 = it1; it2 != factorKeys.end(); it2++) {
        //         Key higherKey = *it2;
        //         // New factors are not part of the changedColStructure
        //         // As they have to be loaded regardless of edit or reconstruct
        //         sortedFactorKeys[lowerKey].insert(higherKey);
        //     }
        // }
    }
    jacobian_.resolveAllocate(-1);

    // // Sort node colStructure
    // for(auto& p : sortedFactorKeys) {
    //     Key lowerKey = p.first;
    //     set<Key>& s = p.second;
    //     auto& factorColStructure = nodes_[lowerKey]->factorColStructure;
    //     for(Key& higherKey : factorColStructure) {
    //         assert(!orderingLess(higherKey, lowerKey));
    //         s.insert(higherKey);
    //     }
    //     // Pull the sorted keys out
    //     factorColStructure.clear();
    //     factorColStructure.insert(factorColStructure.end(), s.begin(), s.end());
    // }

    // cout << "Affected keys: " << endl;
    // for(Key k : *affectedKeys) {
    //     cout << k << " ";
    // }
    // cout << endl;

    // checkInvariant_afterMarkAffected();
}

// Mark all ancestors of directly changed keys, disconnect child from parent and from clique
void CholeskyEliminationTree::markAncestors(const KeySet& affectedKeys, KeySet* markedKeys) {
    // cout << "[CholeskyEliminationTree] markAncestors()" << endl;
    for(const Key key : affectedKeys) {
        markKey(key, markedKeys);
    }
    for(Key k : *markedKeys) {  // marked keys are always backsolved
        backSolveKeys_[k] = true;
    }

    // cout << "Marked keys: ";
    // for(Key k : *markedKeys) {
    //     cout << k << " ";
    // }
    // cout << endl;

    // checkInvariant_afterMarkAncestors();

}

void CholeskyEliminationTree::markKey(const Key key, KeySet* markedKeys) {
    if(nodes_[key]->clique->marked) {
        // Node is already processed
        return;
    }

    // cout << "[CholeskyEliminationTree] markKey() " << key << endl;
    sharedNode node = nodes_[key];

    sharedClique curClique = node->clique;
    Key curKey = key;
    do {
        curClique = curClique->markClique(curKey, markedKeys);
        if(curClique) {
            curKey = curClique->frontKey();
        }
    } while(curClique != nullptr);
}

size_t CholeskyEliminationTree::colWidth(const Key key) const {
    return cholesky_.at(key).width();
}

void CholeskyEliminationTree::symbolicElimination(const KeySet& markedKeys) {
    // cout << "[CholeskyEliminationTree] symbolicElimination()" << endl;

    // TODO: Just add the marked keys in sorted order
    vector<Key> sortedMarkedKeys;
    sortedMarkedKeys.reserve(markedKeys.size());
    sortedMarkedKeys.insert(sortedMarkedKeys.begin(), markedKeys.begin(), markedKeys.end());
    sort(sortedMarkedKeys.begin(), sortedMarkedKeys.end(), orderingLess);
    for(const Key key : sortedMarkedKeys) {
        symbolicEliminateKey(key);
    }

    // // DEBUG Checks
    // for(sharedNode node : nodes_) {
    //     const auto& column = cholesky_.column(node->key);
    //     if(!column.allocated()) {
    //         cout << "After symblic Column " << node->key << " not allocated!" << endl;
    //         exit(1);
    //     }
    // }

    // checkInvariant_afterSymbolic();
}

void CholeskyEliminationTree::symbolicEliminateKey(const Key key) {
    // cout << "[CholeskyEliminationTree] symbolicEliminateKey: " << key << endl;
    sharedNode node = nodes_[key];
    sharedClique clique = node->clique;
    // cout << "clique: ";
    // for(sharedNode node : clique->nodes) {
    //     cout << node->key << " ";
    // }
    // cout << endl;

    if(clique->orderingVersion() != ordering_version_) {
        reorderClique(clique);
    }

    assert(clique->nodes.size() == 1);
    assert(clique->marked == true);

    // Add keys induced by raw factors but only keys that are higher than this key
    node->colStructure.clear();
    node->colStructure.insert(node->colStructure.end(),
                              node->factorColStructure.find(node->key),
                              node->factorColStructure.end());
    assert(sorted_no_duplicates(node->colStructure));

    // cout << " factor col structure: ";
    // for(Key k : node->colStructure) {
    //     cout << k << " ";
    // }

    for(sharedClique childClique : clique->children) {
        if(childClique->orderingVersion() != ordering_version_) {
            reorderClique(childClique);
        }
        sharedNode child = childClique->back();
        assert(sorted_no_duplicates(child->colStructure));
        mergeChildColStructure(node, child);
        assert(sorted_no_duplicates(node->colStructure));
    }

    // cout << "col structure: ";
    // for(Key k : node->colStructure) {
    //     cout << k << " ";
    // }
    // cout << endl;

    // Find parent
    clique->findParent();

    // Check supernode. Only merge if we only have one child and child is marked
    if(clique->children.size() == 1) {
        sharedClique childClique = *(clique->children.begin());
        sharedNode child = childClique->back();
        assert(node->colStructure.size() >= child->colStructure.size() - 1);
        if(childClique->marked 
                && node->colStructure.size() == child->colStructure.size() - 1) {
            // clique->printClique(cout);
            childClique->mergeClique(clique);
            // Need to update our pointer to the current cique
            assert(node->clique = childClique);
            clique = childClique;
            
        }
    }
    // clique->printClique(cout);

    // Set root after merging
    if(clique->parent == nullptr) {
        assert(ordering_[node->key] == ordering_.size() - 1);
        root_ = clique;
    }

    // Set ancestors' descendants
    for(size_t i = 1; i < node->colStructure.size(); i++) {
        // Don't need to count self
        size_t changeAmount = node->colStructure.size() - i;
        Key& ancestorKey = node->colStructure[i];
        auto& changedDescendants = changedDescendants_.at(ancestorKey);
        assert(changedDescendants.empty()
                // || orderingLess(changedDescendants.back().first, node->key));
                || orderingLess(changedDescendants.back(), node->key));
        // changedDescendants.push_back({node->key, changeAmount});
        changedDescendants.push_back(node->key);
    }

    // for(const Key ancestorKey : node->colStructure) {
    //     auto& changedDescendants = changedDescendants_.at(ancestorKey);
    //     // assert(sorted_no_duplicates(changedDescendants));
    //     assert(changedDescendants.empty()
    //             || orderingLess(changedDescendants.back().first, node->key));
    //     changedDescendants_.at(ancestorKey).push_back({node->key});
    // }

    // Merge changed descendants with descendants
    auto& changedDescendants = changedDescendants_.at(node->key);
    // assert(sorted_no_duplicates(changedDescendants));
    if(!changedDescendants.empty()) {
        size_t i = 0;
        auto& descendants = descendants_.at(node->key);
        // assert(sorted_no_duplicates(descendants));
        for(; i < descendants.size(); i++) {
            // Find the index in descendants that is the same as the 
            // first element in changed descendants. If not found,
            // i will be set to the end of descendants
            if(!orderingLess(descendants[i], changedDescendants[0])) {
            // if(!orderingLess(descendants[i].first, changedDescendants[0].first)) {
                break;
            }
            else {
                // assert(orderingLess(descendants[i].first, changedDescendants[0].first));
                assert(orderingLess(descendants[i], changedDescendants[0]));
            }
        }
        // At this point we just need to insert changedDescendants into descendants
        descendants.resize(i + changedDescendants.size());
        for(size_t j = 0; j < changedDescendants.size(); j++) {
            descendants[i + j] = changedDescendants[j];
        }
        // assert(sorted_no_duplicates(descendants_.at(node->key)));
        // assert(sorted_no_duplicates(changedDescendants_.at(node->key)));
    }
}

void CholeskyEliminationTree::mergeChildColStructure(sharedNode parent, sharedNode child) {
    vector<Key>& col1 = parent->colStructure;
    const vector<Key>& col2 = child->colStructure;
    vector<Key> newColStructure;
    newColStructure.reserve(col1.size() + col2.size());

    size_t i1 = 0, i2 = 1;  // merged col structure does not include diagonal

    while(i1 < col1.size() || i2 < col2.size()) {
        if(i2 == col2.size()) {
            newColStructure.push_back(col1[i1]);
            i1++;
        } 
        else if(i1 == col1.size()) {
            newColStructure.push_back(col2[i2]);
            i2++;
        }
        else if(col1[i1] == col2[i2]) {
            newColStructure.push_back(col1[i1]);
            i1++;
            i2++;
        }
        else if(orderingLess(col1[i1], col2[i2])) {
            newColStructure.push_back(col1[i1]);
            i1++;
        }
        else if(orderingLess(col2[i2], col1[i1])) {
            newColStructure.push_back(col2[i2]);
            i2++;
        }
        else {
            assert(0);
        }
    }
    parent->colStructure = std::move(newColStructure);
}

void CholeskyEliminationTree::constructCSCMatrix(
        const vector<Key>& reorderKeys,
        const KeySet& observedKeys,
        int* nEntries,
        int* nVars,
        int* nFactors,
        vector<Key>* keyOrdering,
        vector<int>* A,
        vector<int>* p,
        vector<int>* cmember) {

    set<FactorIndex> rawFactors;
    set<Key> fixedKeys;
    vector<vector<FactorIndex>> cscMatrix(reorderKeys.size(), vector<FactorIndex>());
    unordered_map<Key, size_t> keyMap;

    int keyCount = 0;
    for(const Key key : reorderKeys) {
        keyMap[key] = keyCount;
        keyCount++;
        sharedNode node = nodes_[key];
        for(FactorIndex factorIndex : node->factorIndices) {
            rawFactors.insert(factorIndex);
        }
        for(sharedClique childClique : node->clique->children) {
            assert(!childClique->marked);
            assert(is_reordered_[childClique->backKey()] == false);
            assert(childClique->back()->colStructure[1] == node->key);
            fixedKeys.insert(childClique->backKey());

        }
    }

    size_t factorCount = 0;
    for(FactorIndex factorIndex : rawFactors) {
        sharedFactor factor = factors_[factorIndex]; 
        for(const Key key : factor->keys()) {
            if(is_reordered_[key]) {
                cscMatrix.at(keyMap.at(key)).push_back(factorCount);
            }
        }
        factorCount++;
    }

    for(Key key : fixedKeys) {
        sharedNode node = nodes_[key];
        for(size_t i = 1; i < node->colStructure.size(); i++) {
            Key colKey = node->colStructure[i];

            assert(!node->clique->marked);
            assert(node->clique->parent->marked);
            assert(node->colStructure[1] == node->clique->parent->frontKey());
            assert(is_reordered_[node->clique->parent->frontKey()]);
            assert(is_reordered_[colKey] == true);

            cscMatrix.at(keyMap.at(colKey)).push_back(factorCount);
        }
        factorCount++;

    }

    p->reserve(reorderKeys.size() + 1);
    p->push_back(0);

    size_t observedKeyConstraint = (reorderKeys.size() == observedKeys.size())? 0 : 1;
    size_t count = 0;

    for(size_t i = 0; i < reorderKeys.size(); i++) {
        Key key = reorderKeys[i];
        if(observedKeys.find(key) != observedKeys.end()) {
            cmember->push_back(observedKeyConstraint);
        }
        else {
            cmember->push_back(0);
        }
        for(FactorIndex factorIndex : cscMatrix[i]) {
            A->push_back(factorIndex);
            count++;
        }
        p->push_back(count);
    }

    *nEntries = count;
    *nFactors = factorCount;
    *nVars = reorderKeys.size();

}

/*
void CholeskyEliminationTree::constructCSCMatrix(
        const vector<Key>& reorderKeys,
        const KeySet& observedKeys,
        int* nEntries,
        int* nVars,
        int* nFactors,
        vector<Key>* keyOrdering,
        vector<int>* A,
        vector<int>* p,
        vector<int>* cmember) {

    map<Key, set<FactorIndex>> factorMap;
    unordered_map<FactorIndex, int> allFactors;
    set<Key, OrderingLess> fixedKeys(orderingLess);

    for(const Key key : reorderKeys) {
        sharedNode node = nodes_[key];
        for(FactorIndex factorIndex : node->factorIndices) {
            // map factor to an index that starts from 0
            allFactors.insert({factorIndex, allFactors.size()});

            sharedFactor factor = factors_[factorIndex]; 
            for(const Key otherKey : factor->keys()) {
                factorMap[otherKey].insert(allFactors.at(factorIndex));
                if(!is_reordered_[otherKey]) {
                    fixedKeys.insert(otherKey);
                }
            }
        }
    }

    int count = 0;

    // cout << "fixed keys size = " << fixedKeys.size() << endl;

    p->reserve(factorMap.size() + 1);
    p->push_back(0);
    // Do fixedKeys first
    for(const Key key : fixedKeys) {
        keyOrdering->push_back(key);
        cmember->push_back(cmember->size());
        const set<FactorIndex>& s = factorMap.at(key);
        for(const FactorIndex factorIndex : s) {
            A->push_back(factorIndex);
            count++;
        }
        p->push_back(count);
    }

    size_t constraintKey = cmember->size();
    size_t observedConstraintKey = (reorderKeys.size() == observedKeys.size())? 
                                        constraintKey : constraintKey + 1;

    // Then do reorderKeys 
    for(const Key key : reorderKeys) {
        keyOrdering->push_back(key);
        if(observedKeys.find(key) != observedKeys.end()) {
            cmember->push_back(observedConstraintKey);
        }
        else {
            cmember->push_back(constraintKey);
        }
        const set<FactorIndex>& s = factorMap.at(key);
        for(const FactorIndex factorIndex : s) {
            A->push_back(factorIndex);
            count++;
        }
        p->push_back(count);
    }


    *nEntries = count;
    *nFactors = allFactors.size();
    *nVars = factorMap.size();

}
*/

/*
void CholeskyEliminationTree::constructCSCMatrix(
        const vector<Key>& reorderKeys,
        int* nEntries,
        int* nVars,
        int* nFactors,
        vector<int>* A,
        vector<int>* p,
        vector<int>* cmember,
        unordered_set<Key>* is_reordered) {
    map<Key, set<FactorIndex>> factorMap;
    unordered_map<FactorIndex, int> allFactors;
    unordered_set<Key> fixedKeys;

    for(const Key key : reorderKeys) {
        is_reordered->insert(key);
    }

    for(const Key key : reorderKeys) {
        sharedNode node = nodes_[key];
        for(FactorIndex factorIndex : node->factorIndices) {
            // map factor to an index that starts from 0
            allFactors.insert({factorIndex, allFactors.size()});

            sharedFactor factor = factors_[factorIndex]; 
            for(const Key otherKey : factor->keys()) {
                factorMap[otherKey].insert(factorIndex);
                if(is_reordered->find(otherKey) == is_reordered->end()) {
                    fixedKeys.insert(otherKey);
                }
            }
        }
    }

    int count = 0;

    p->reserve(factorMap.size() + 1);
    p->push_back(0);
    // Do reorderKeys first
    for(const Key key : reorderKeys) {
        const set<FactorIndex>& s = factorMap.at(key);
        for(const FactorIndex factorIndex : s) {
            A->push_back(allFactors[factorIndex]);
            count++;
        }
        p->push_back(count);
    }


    // Then do fixed keys with constraint 0
    for(const Key key : fixedKeys) {
        const set<FactorIndex>& s = factorMap.at(key);
        for(const FactorIndex factorIndex : s) {
            A->push_back(allFactors[factorIndex]);
            count++;
        }
        p->push_back(count);
    }

    // Use constraints to fix order of fixed keys against reorder keys
    // fixed keys should be constrained to be in front of reorder keys
    cmember->resize(factorMap.size());

    // Set up fixed keys constraints
    // Need to do this kind of convoluted thing to make sure there are no gaps in cmember
    size_t constraintKey = 0, newConstraintKey = 0;
    for(int i = 0; i < fixedKeys.size(); i++) {
        cmember->at(cmember->size() - 1 - i) = constraintKey;
        newConstraintKey = constraintKey + 1;
    }
    constraintKey = newConstraintKey;

    // Set up reorder keys constraints
    for(int i = 0; i < reorderKeys.size() - 2; i++) {
        cmember->at(i) = constraintKey;
        newConstraintKey = constraintKey + 1;
    }
    constraintKey = newConstraintKey;

    // Set up the last two keys constraints
    for(int i = reorderKeys.size() - 2; i < reorderKeys.size() - 1; i++) {
        cmember->at(i) = constraintKey;
        newConstraintKey = constraintKey + 1;
    }
    constraintKey = newConstraintKey;
    for(int i = reorderKeys.size() - 1; i < reorderKeys.size(); i++) {
        cmember->at(i) = constraintKey;
        newConstraintKey = constraintKey + 1;
    }


    *nEntries = count;
    *nFactors = allFactors.size();
    *nVars = factorMap.size();
}
*/

void CholeskyEliminationTree::getPartialReordering(const vector<Key>& reorderKeys,
                                                   const KeySet& observedKeys,
                                                   vector<Key>* partialOrdering) {
    if(reorderKeys.size() == 0) {
        return;
    }
    else if(reorderKeys.size() == 1) {
        partialOrdering->push_back(reorderKeys.front());
        return;
    }

    int nEntries, nVars, nFactors;
    vector<int> A, p, cmember;
    vector<Key> keyOrdering;
    // unordered_set<Key> is_reordered;

    constructCSCMatrix(reorderKeys, observedKeys,
                       &nEntries, &nVars, &nFactors, 
                       &keyOrdering, &A, &p, &cmember);

    assert(p.size() == nVars + 1);

    const size_t Alen = ccolamd_recommended((int) nEntries, (int) nFactors, (int) nVars);

    A.resize(Alen);

    //double* knobs = nullptr; /* colamd arg 6: parameters (uses defaults if nullptr) */
    double knobs[CCOLAMD_KNOBS];
    ccolamd_set_defaults(knobs);
    knobs[CCOLAMD_DENSE_ROW] = -1;
    knobs[CCOLAMD_DENSE_COL] = -1;

    int stats[CCOLAMD_STATS]; /* colamd arg 7: colamd output statistics and error codes */

    // cout << "A: ";
    // for(int i : A) {
    //     cout << i << " ";
    // }
    // cout << endl;
    // cout << "p: ";
    // for(int i : p) {
    //     cout << i << " ";
    // }
    // cout << endl;
    // cout << "cmember: ";
    // for(int i : cmember) {
    //     cout << i << " ";
    // }
    // cout << endl;
    

    // cout << "p: ";
    // for(int i = 0; i < p.size(); i++) {
    //     cout << p.at(i) << " ";
    // }
    // cout << endl << endl;

    // call colamd, result will be in p
    /* returns (1) if successful, (0) otherwise*/
    if (nVars > 0) {
        gttic(ccolamd);
        int rv = ccolamd((int) nFactors, (int) nVars, (int) Alen, &A[0], &p[0],
                knobs, stats, &cmember[0]);
        if (rv != 1) {
            throw runtime_error("ccolamd failed with return value " + to_string(rv));
        }
    }

    // cout << "p: ";
    // for(int i = 0; i < p.size(); i++) {
    //     cout << p.at(i) << " ";
    // }
    // cout << endl << endl;

    partialOrdering->resize(reorderKeys.size());
    size_t index = 0;
    // cout << "p: ";
    // for(int i = 0; i < p.size(); i++) {
    //     cout << p[i] << " ";
    // }
    // cout << endl;
    for(int i = 0; i < p.size(); i++) {
        if(p[i] == -1) {
            break;
        }
        Key key = reorderKeys[p[i]];
        if(is_reordered_[key]) {
            partialOrdering->at(index++) = key;
        }
        // if(p[i] >= reorderKeys.size()) {
        //     continue;
        // }
        // else if(p[i] == -1) {
        //     break;
        // }
        // const Key key = reorderKeys[p[i]];
        // // if(is_reordered.find(key) != is_reordered.end()) {
        //     // We only care about the keys in reorderKeys
        //     partialOrdering->at(index++) = key;
        // // }
    }

    // if(nodes_.size() < 300) {
    //     Ordering ordering(*partialOrdering);
    //     cout << "Ordering: ";
    //     ordering.print();
    //     cout << endl;
    // }
}

void CholeskyEliminationTree::updateOrdering(const KeySet& markedKeys, 
                                             const KeySet& observedKeys) {
    // cout << "[CholeskyEliminationTree] updateOrdering()" << endl;
    vector<Key> reorderKeys;
    vector<Key> partialOrdering;

    reorderKeys.insert(reorderKeys.begin(), markedKeys.begin(), markedKeys.end());

    for(Key key : reorderKeys) {
        assert(nodes_[key]->clique->marked);
        assert(is_reordered_[key] == false);
        is_reordered_[key] = true;
    }

    getPartialReordering(reorderKeys, observedKeys, &partialOrdering);

    // Adjust ordering_. Shift all fixed keys to the front
    size_t lowestReorderedIndex = -1;
    vector<Key> deltaReorderKeys;
    vector<Key> newOrderingToKey;
    newOrderingToKey.reserve(orderingToKey_.size());
    for(size_t i = 0; i < ordering_.size(); i++) {
        Key key = orderingToKey_[i];
        if(!is_reordered_[key]) {
            newOrderingToKey.push_back(key);
            if(lowestReorderedIndex != -1) {
                deltaReorderKeys.push_back(key);       
            }
        }
        else if(lowestReorderedIndex == -1) {
            lowestReorderedIndex = i;
        }
    }
    newOrderingToKey.insert(newOrderingToKey.end(), 
                            partialOrdering.begin(), 
                            partialOrdering.end());

    assert(newOrderingToKey.size() == ordering_.size());

    deltaReorderKeys.insert(deltaReorderKeys.end(), 
                            partialOrdering.begin(), 
                            partialOrdering.end());

    // cout << "Delta reorder keys: ";
    // for(Key k : deltaReorderKeys) {
    //     cout << k << " ";
    // }
    // cout << endl;
    // cout << "Partial ordering: ";
    // for(Key k : partialOrdering) {
    //     cout << k << " ";
    // }
    // cout << endl;
    
    for(int i = 0; i < newOrderingToKey.size(); i++) {
        Key key = newOrderingToKey[i];
        ordering_[key] = i;
    }

    orderingToKey_ = std::move(newOrderingToKey);

    // size_t curIndex = 0;
    // for(size_t i = 0; i < orderingToKey_.size(); i++) {
    //     Key key = orderingToKey_[i];
    //     if(markedKeys.find(key) == markedKeys.end()) {

    //         if(lowestReorderedIndex == -1 && ordering_[key] != curIndex) {
    //             lowestReorderedIndex = curIndex;
    //             // cout << "lowestReorderedIndex = " << lowestReorderedIndex << endl;
    //         }

    //         orderingToKey_[curIndex] = key;
    //         ordering_[key] = curIndex;
    //         // cout << "Key " << key << " is now position " << curIndex << endl;;

    //         curIndex++;
    //     }
    //     else {
    //         if(lowestReorderedIndex == -1) {
    //             lowestReorderedIndex = i;
    //             // cout << "lowestReorderedIndex = " << lowestReorderedIndex << endl;
    //         }
    //     }
    // }

    ordering_version_++;    // update ordering version

    if(!delta_.allocated()) {
        cout << "delta not allocated!" << endl;
        exit(1);
    }
    // cout << "reorder keys: ";
    // for(Key k : partialOrdering) {
    //     cout << k << " ";
    // }
    // cout << endl;
    delta_.reorderBlocks(deltaReorderKeys, lowestReorderedIndex);


    // cout << "update ordering here4" << endl;

    // for(size_t i = 0; i < partialOrdering.size(); i++) {
    //     Key key = partialOrdering[i];
    //     orderingToKey_[curIndex + i] = key;
    //     ordering_[key] = curIndex + i;
    //     // cout << "Key " << key << " is now position " << curIndex + i << endl;
    // }

    // Reset all affected nodes
    for(size_t i = 0; i < partialOrdering.size(); i++) {
        Key key = partialOrdering[i];
        sharedNode node = nodes_[key];
        node->colStructure.clear();

        // If the node has any children cliques, the children need to be reattached
        // Since they might have new parents
        // Need to copy as changing node->clique->children invalidate iterators
        vector<sharedClique> cliqueChildren;
        for(sharedClique child : node->clique->children) {
            cliqueChildren.push_back(child);
        }
        for(sharedClique child : cliqueChildren) {
            reparentOrphanClique(child);
        }

        // node->colContribution.clear();
        // assert(node->changedColStructure.empty());
        // assert(node->changedColContribution.empty());
        // assert(node->newColContribution.empty());

        assert(node->ordering_version != ordering_version_);
        node->ordering_version = ordering_version_;

        // We had a problem here with cliques not inheriting reconstructCols correctly
        // after reset, since inheritCols requires a sorted version of the 
        // old keys in the new ordering. For edits and for reconstruct (no reordering), 
        // this is not a problem because there is guaranteed to be no reordering, and 
        // blockStartVec provides the old column structure
        // For reconstructs (after reordering) we need the old keys
        cholesky_.resetColumn(key);
        // cout << "reset column " << key << endl;

        node->is_reordered = true;

        // Need to sort factor col structure
        set<Key, OrderingLess> newFactorColStructure(orderingLess);
        set<Key, OrderingLess> newChangedFactorColStructure(orderingLess);
        for(Key k : node->factorColStructure) {
            newFactorColStructure.insert(k);
        }
        node->factorColStructure = std::move(newFactorColStructure);
        for(Key k : node->changedFactorColStructure) {
            if(!orderingLess(k, node->key)) {
                newChangedFactorColStructure.insert(k);
            }
        }
        node->changedFactorColStructure = std::move(newChangedFactorColStructure);


        // Need to also change descendants 
        /*
        node->factorColStructure.clear();
        node->changedFactorColStructure.clear();
        for(const FactorIndex factorIndex : node->factorIndices) {
            sharedFactor factor = factors_[factorIndex];
            for(const Key otherKey : factor->keys()) {
                if(!orderingLess(otherKey, key)) {
                    // otherKey is greater than or equal to us
                    node->factorColStructure.push_back(otherKey);
                    if(factorLinearizeStatus_[factorIndex] == RELINEARIZED) {
                        node->changedFactorColStructure.push_back(otherKey);
                    }
                }
            }
        }

        std::sort(node->factorColStructure.begin(), 
                  node->factorColStructure.end(),
                  orderingLess);

        std::sort(node->changedFactorColStructure.begin(), 
                  node->changedFactorColStructure.end(),
                  orderingLess);
        */

        // Fix descendants
        // Descendants in the reordered keys need to be removed
        // While descendants that are not marked/reordered need to be kept
        // since we're not running symbolic elim on them
        assert(changedDescendants_.at(key).empty());
        auto& descendants = descendants_.at(key);
        for(size_t j = 0; j < descendants.size(); j++) {
            if(is_reordered_[descendants[j]]) {
            // if(is_reordered_[descendants[j].first]) {
                descendants.resize(j);  // Resize everything after the first reordered key
                break;
            }
        }
    }
}

void CholeskyEliminationTree::reparentOrphanClique(sharedClique clique) {
    // cout << "in reparent orphan" << endl;
    reorderClique(clique);

    clique->findParent();

    assert(clique->parent != nullptr);
}

/*
void CholeskyEliminationTree::updateOrdering(KeySet* affectedKeys) {
    // cout << "[CholeskyEliminationTree] updateOrdering()" << endl;
    getTotalReordering();

    // TODO: Support partial ordering

    // For each variable, reset all relevant data structures and check all factors
    y_.resetBlocks(false);
    delta_.resetBlocks(false);

    map<size_t, Key> orderingToKey_;
    for(size_t key = 0; key < ordering_.size(); key++) {
        orderingToKey_.insert({ordering_[key], key});
        sharedNode node = nodes_[key];
        node->colStructure.clear();
        node->colContribution.clear();
        assert(node->changedColStructure.empty());
        assert(node->changedColContribution.empty());
        assert(node->newColContribution.empty());
        assert(node->relinearize == false);
        assert(node->marked == false);

        affectedKeys->insert(key);
        node->marked = true;

        cholesky_.resetColumn(key);

        // Re-add all factors
        for(const FactorIndex factorIndex : node->factorIndices) {
            assert(factorLinearizeStatus_[factorIndex] == LINEARIZED);
            sharedFactor factor = factors_[factorIndex];
            for(const Key otherKey : factor->keys()) {
                if(!orderingLess(otherKey, key)) {
                    // If the other key is in our column, it adds to a contribution block
                    // and a column structure block
                    node->newColContribution[key]++;
                    node->changedColStructure.insert({otherKey, 1});
                }
                
            }
        }
    }


    for(int i = 0; i < ordering_.size(); i++) {
        const Key key = orderingToKey_[i];

        size_t width = cholesky_.column(key).width();
        y_.preallocateBlock(key, width, true);
        delta_.preallocateBlock(key, width, true);
    }

    y_.resolveAllocate();
    delta_.resolveAllocate();
} 
*/

void CholeskyEliminationTree::setEditOrReconstruct(sharedClique clique) {

    // If reordered, automatically reconstruct
    if(clique->front()->is_reordered) {
    // if(true) {
        clique->is_reconstruct = true;
    }
    else {
        size_t totalCols = 0;
        size_t changedCols = 0;

        for(sharedNode node : clique->nodes) {
            Key& nodeKey = node->key;
            auto it1 = node->factorColStructure.find(nodeKey);
            while(it1 != node->factorColStructure.end()) {
                totalCols++;
                it1++;
            }
            auto it2 = node->changedFactorColStructure.find(nodeKey);
            while(it2 != node->changedFactorColStructure.end()) {
                changedCols++;
                it2++;
            }
            // for(auto& p : descendants_.at(nodeKey)) {
            //     totalCols += p.second;
            // }
            // for(auto& p : changedDescendants_.at(nodeKey)) {
            //     changedCols += p.second;
            // }
            totalCols += descendants_.at(nodeKey).size();
            changedCols += changedDescendants_.at(nodeKey).size();
            // for(const Key key : descendants_.at(nodeKey)) {
            //     cliqueDescendants.insert(key);
            // }
            // for(const Key key : changedDescendants_.at(nodeKey)) {

            //     cliqueChangedDescendants.insert(key);
            // }
        }

        if(totalCols == 0 || changedCols >= 0.9 * totalCols) {
            clique->is_reconstruct = true;
        }
        else {
            clique->is_reconstruct = false;
        }
    }

}

/*
void CholeskyEliminationTree::setEditOrReconstruct(sharedClique clique) {

    // If reordered, automatically reconstruct
    if(clique->front()->is_reordered) {
    // if(true) {
        clique->is_reconstruct = true;
    }
    else {

        unordered_set<Key> cliqueFactorColStructure;
        unordered_set<Key> cliqueChangedFactorColStructure;
        unordered_set<Key> cliqueDescendants;
        unordered_set<Key> cliqueChangedDescendants;
        // TODO: optimize for one node cliques
        for(sharedNode node : clique->nodes) {
            Key& nodeKey = node->key;
            auto it1 = node->factorColStructure.find(nodeKey);
            while(it1 != node->factorColStructure.end()) {
                cliqueFactorColStructure.insert(*it1);
                it1++;
            }
            auto it2 = node->changedFactorColStructure.find(nodeKey);
            while(it2 != node->changedFactorColStructure.end()) {
                cliqueChangedFactorColStructure.insert(*it2);
                it2++;
            }
            // for(const Key key : node->factorColStructure) {
            //     cliqueFactorColStructure.insert(key);
            // }
            // for(const Key key : node->changedFactorColStructure) {
            //     cliqueChangedFactorColStructure.insert(key);
            // }
            for(const Key key : descendants_.at(nodeKey)) {
                cliqueDescendants.insert(key);
            }
            for(const Key key : changedDescendants_.at(nodeKey)) {
                cliqueChangedDescendants.insert(key);
            }
        }

        // Remove clique keys
        for(sharedNode node : clique->nodes) {
            const Key key = node->key;
            cliqueFactorColStructure.erase(key);
            cliqueChangedFactorColStructure.erase(key);
            cliqueDescendants.erase(key);
            cliqueChangedDescendants.erase(key);
        }

        size_t totalCols = cliqueFactorColStructure.size() + cliqueDescendants.size();
        size_t changedCols = cliqueChangedFactorColStructure.size() 
            + cliqueChangedDescendants.size();
        if(totalCols == 0 || changedCols >= 0.95 * totalCols) {
            clique->is_reconstruct = true;
        }
        else {
            clique->is_reconstruct = false;
        }
    }

}
}
*/

void CholeskyEliminationTree::choleskyElimination(const Values& theta) {
    // cout << "[CholeskyEliminationTree] choleskyElimination()" << endl;
    vector<pair<sharedClique, bool>> stack(1, {root_, false});
    while(!stack.empty()) {
        auto& curPair = stack.back();
        sharedClique clique = curPair.first;
        sharedClique parent = clique->parent;
        bool& expanded = curPair.second;
        if(!expanded) {
            expanded = true;

            // cout << "Restore pass: ";
            // for(sharedNode node : clique->nodes) {
            //     cout << node->key << " ";
            // }
            // cout << endl;

            // Need to reorder fixed nodes before inheriting
            bool has_reconstruct = false;
            if(!clique->marked) {
                if(clique->orderingVersion() != ordering_version_) {
                    // Reorder each node in the clique
                    reorderClique(clique);
                }
            }

            // Only expand to children if we are marked or have a reconstruct column
            // It is possible to be marked but all ancestors are edits
            if(clique->marked) {
                assert(clique->orderingVersion() == ordering_version_);

                // Decide if edit or reconstruct
                setEditOrReconstruct(clique);

                // edit from this clique and restore to linear state if edit
                editAndRestoreFromClique(clique);


                if(clique->is_reconstruct) {
                    // Need to add to reconstructCols and editCols in reverse order
                    for(sharedNode node : clique->nodes) {
                        if(markedStatus_[node->key] != NEW) {
                            cholesky_.resetColumn(node->key);
                            assert(markedStatus_[node->key] == UNMARKED);
                            markedStatus_[node->key] = RECONSTRUCT;
                            // cout << "Set " << node->key << " to RECONSTRUCT" << endl;
                        }
                    }
                }
                else {
                    for(sharedNode node : clique->nodes) {
                        markedStatus_[node->key] = EDIT;
                        // cout << "Set " << node->key << " to EDIT" << endl;
                    }
                    // subtract old linearized blocks
                    prepareEditClique(clique);
                }

                // Allocate new blocks for edit or reconstruct
                for(sharedNode node : clique->nodes) {
                    const Key key = node->key;
                    vector<pair<Key, size_t>> colStructureHeight;
                    colStructureHeight.reserve(node->colStructure.size());
                    for(const auto& otherKey : node->colStructure) {
                        assert(!orderingLess(otherKey, key));
                        colStructureHeight.push_back({otherKey, colWidth(otherKey)});
                    }
                    cholesky_.column(key).preallocateBlocks(colStructureHeight);
                    jacobian_.resolveAllocate(key);
                    cholesky_.resolveAllocate(key);
                }
            }
            else {
                assert(clique->orderingVersion() == ordering_version_);
                has_reconstruct = reconstructFromClique(clique);
            }
            // only expand to children if we are marked or have a reconstruct
            // It is possible to be marked but all ancestors are edits
            if(clique->marked || has_reconstruct) {
                for(sharedClique child : clique->children) {
                    stack.push_back({child, false});
                }
            }

        }
        else {

            // cout << "Eliminate pass: ";
            // for(sharedNode node : clique->nodes) {
            //     cout << node->key << " ";
            // }
            // cout << endl;
            // cout << "is reconstruct = " << clique->is_reconstruct << endl;

            // Eliminate pass
            stack.pop_back();

            // Unmarked clique already reconstructed
            if(clique->marked) {
                // All factors of this clique should be relinearized
                // Do AtA for each node
                prepareEliminateClique(clique, theta);

                // Eliminiate clique
                eliminateClique(clique);

                // Reset node member variables
                for(sharedNode node : clique->nodes) {
                    node->relinearize = false;
                    node->is_reordered = false;
                    node->changedFactorColStructure.clear();
                    changedDescendants_[node->key].clear();
                    markedStatus_[node->key] = UNMARKED;
                    is_reordered_[node->key] = false;
                }
            }

            // Reset member variables
            clique->marked = false;
            clique->is_reconstruct = true;
        }
    }

    // checkInvariant_afterCholesky();

    // // DEBUG Checks
    // for(sharedNode node : nodes_) {
    //     const auto& column = cholesky_.column(node->key);
    //     if(!column.allocated()) {
    //         cout << "Column " << node->key << " not allocated!" << endl;
    //         exit(1);
    //     }
    // }
}

/*
void CholeskyEliminationTree::choleskyElimination(Values& theta) {
    // cout << "[CholeskyEliminationTree] choleskyElimination()" << endl;
    vector<pair<sharedNode, bool>> stack(1, {root_, false});
    while(!stack.empty()) {
        auto& curPair = stack.back();
        sharedNode node = curPair.first;
        sharedNode parent = node->parent;
        const Key key = node->key;
        bool& expanded = curPair.second;
        if(!expanded) {
            // restore pass
            // After restore pass, the node should:
            //    1) Be done with all necessary edits
            //    2) Decide if reconstruct or edit
            //    2) Memory should be allocated
            //    3) If reconstruct, columns should be set to zero and blocks reordered
            //    4) If edit, diagonal block should be multiplied to the column
            expanded = true;
            // cout << "Restore pass key: " << key << endl;

            // Pass down the reconstruct columns and edit columns (for marked nodes)
            assert(node->reconstructCols.empty());
            if(parent) {
                for(const Key ancestorKey : parent->reconstructCols) {
                    if(node->colStructure.find(ancestorKey) != node->colStructure.end()) {
                        assert(orderingLess(key, ancestorKey));
                        node->reconstructCols.push_back(ancestorKey);
                        // cout << "Reconstruct ancestor " << ancestorKey << endl;
                    }
                }

                if(node->marked) {
                    for(const Key ancestorKey : parent->editCols) {
                        assert(orderingLess(key, ancestorKey));
                        if(node->colStructure.find(ancestorKey) != node->colStructure.end()) {
                            node->editCols.push_back(ancestorKey);
                        }
                    }
                }
            }

            // Only expand to children if we are marked or have a reconstruct
            // It is possible to be marked but all ancestors are edits
            if(node->marked || !node->reconstructCols.empty()) {
                if(node->marked) {
                    // FIXME: this is really slow, let's fix
                    handleEdits(node);
                
                    if(node->is_reconstruct) {
                        node->reconstructCols.push_back(key);

                        cholesky_.resetBlocks(node->key);

                    }
                    else {
                        // cout << "node marked. push to edit col" << endl;
                        node->editCols.push_back(key);
                        node->is_reconstruct = false;

                        // need to restore column and subtract the old linearized blocks
                        prepareEditColumn(node);
                    }
                    
                    // Allocate blocks for edit or reconstruct
                    for(const auto& p : node->colStructure) {
                        const Key otherKey = p.first;
                        assert(!orderingLess(otherKey, key));
                        // If reconstruct, initialize blocks
                        bool alloc = cholesky_.preallocateBlock(otherKey, key, node->is_reconstruct);
                    }

                    jacobian_.resolveAllocate(node->key);
                    cholesky_.resolveAllocate(node->key);
                }
                else {
                    // handle reconstructs
                    handleReconstructs(node);
                }

                // Only expand to children if we are marked or have a reconstruct
                // It is possible to be marked but all ancestors are edits
                for(sharedNode child : node->children) {
                    stack.push_back({child, false});
                }

            }

            // Set edit or reconstruct
            // Compare edit cost and reconstruct cost
            // Edit involves (extra cost):
            //  1) Multiplying diagonal block to the column (len(L) mults)
            //  2) Restoring the blocks that were changed (len(changeContribution) adds + mults)
            //  Everything else remains the same
            // Reconstruct involves (extra cost):
            //  1) Setting all to 0 (memset; doesn't count)
            //  2) Reconstructing all the blocks that were not changed but were reset (len(notChangedContribution) adds + mults)
            //  So roughly, if 2*len(changedContribution) < len(totalContribution) - len(L)
            //  We should do edit
        }
        else {
            // cout << "Eliminate pass key: " << key << endl;
            // eliminate pass
            // Appropriate data structure should be reset
            stack.pop_back();

            // Unmarked nodes already reconstructed
            if(node->marked) {
                // cout << "AtA" << endl;

                // All factors of this node should be relinearized
                for(const FactorIndex factorIndex : node->factorIndices) {
                    sharedFactor factor = factors_[factorIndex];
                    if(!node->is_reconstruct && factorLinearizeStatus_[factorIndex] == LINEARIZED) {
                        // Can skip factor if node is edit and factor has already been linearized
                        continue;
                    }
                    if(factorLinearizeStatus_[factorIndex] != LINEARIZED) {
                        vector<Matrix> A;
                        Vector b(factor->dim());
                        factor->linearizeToMatrix(theta, &A, &b);
                        for(int i = 0; i < factor->size(); i++) {
                            const Key factorKey = factor->keys()[i];
                            jacobian_.block(factorIndex, factorKey) = A[i];
                            b_.block(factorIndex) = b;
                        }
                        factorLinearizeStatus_[factorIndex] = LINEARIZED;
                    }
                    // Do AtA and Atb here, but just for this column
                    // If edit, only need to do factors that are un/relienarized
                    SparseColumnBlockMatrix& column = cholesky_.column(key);
                    for(const Key otherKey : factor->keys()) {
                        if(!orderingLess(otherKey, key)) {
                            assert(column.blockExists(otherKey));
                            // TODO: AtA can probably be done at once for each row
                            column.block(otherKey) += jacobian_.block(factorIndex, otherKey).transpose() * jacobian_.block(factorIndex, key);
                        }
                    }
                    y_.block(key) += jacobian_.block(factorIndex, key).transpose() * b_.block(factorIndex);
                }

                // If in a supernode and adopted, then can skip this step. Ancestor will eliminate
                // If in a supernode but is at the top of supernode, construct the column to eliminate
                // If not in a supernode, then just eliminate
                if(!node->adoptedSkip) {
                    // If not in supernode or if at top of supernode
                    if(node->adoptedCols.empty()) {
                        // if not in supernode
                        // Do elimination on this node
                        const size_t colWidth = cholesky_.column(key).width();
                        Eigen::VectorXd diagonalY(colWidth);
                        // copy data over
                        diagonalY = y_.block(key);
                        eliminateColumn(&cholesky_.column(key), &diagonalY);
                        
                        // copy data back
                        y_.block(key) = diagonalY;
                    }
                    else {
                        // construct supernode

                        // technically this node also adopts self
                        node->adoptedCols.push_back(key);

                        // We only need to allocate super delta for diagonal
                        // Aggregate all widths
                        size_t totalWidth = 0;
                        for(const Key adoptedKey : node->adoptedCols) {
                            const size_t adoptedWidth = cholesky_.column(adoptedKey).width();
                            totalWidth += adoptedWidth;
                        }
                        size_t totalHeight = cholesky_.column(node->adoptedCols.front()).height();
                        SparseColumnBlockMatrix superColumn(INT_MAX, totalWidth, true);
                        Eigen::VectorXd diagonalY(totalWidth);
                        auto iter = node->colStructure.begin();
                        iter++; // don't include diagonal as we already have a big diagonal
                        for(; iter != node->colStructure.end(); iter++) {
                            const Key otherKey = iter->first;
                            const Key otherWidth = cholesky_.column(otherKey).width();
                            superColumn.preallocateBlock(otherKey, otherWidth, true);
                        }
                        superColumn.resolveAllocate();
                        assert(superColumn.height() == totalHeight);

                        // Copy data over
                        size_t startCol = 0;
                        for(const Key adoptedKey : node->adoptedCols) {
                            const auto& adoptedColumn = cholesky_.column(adoptedKey);
                            size_t adoptedWidth = adoptedColumn.width();
                            size_t adoptedHeight = adoptedColumn.height();
                            assert(totalHeight - startCol == adoptedHeight);
                            Block colBlock = superColumn.submatrix(startCol, startCol, 
                                                  adoptedHeight, adoptedWidth);
                            colBlock = adoptedColumn.blockRange(0, adoptedHeight);

                            diagonalY.block(startCol, 0, adoptedWidth, 1) = y_.block(adoptedKey);

                            startCol += adoptedWidth;
                        }

                        eliminateColumn(&superColumn, &diagonalY);

                        // Copy data back
                        startCol = 0;
                        for(const Key adoptedKey : node->adoptedCols) {
                            auto& adoptedColumn = cholesky_.column(adoptedKey);
                            size_t adoptedWidth = adoptedColumn.width();
                            size_t adoptedHeight = adoptedColumn.height();
                            assert(totalHeight - startCol == adoptedHeight);
                            Block colBlock = superColumn.submatrix(startCol, startCol, 
                                                  adoptedHeight, adoptedWidth);
                            adoptedColumn.blockRange(0, adoptedHeight) = colBlock;

                            y_.block(adoptedKey) 
                                = diagonalY.block(startCol, 0, adoptedWidth, 1);

                            startCol += adoptedWidth;
                        }

                    }
                }

                // cout << "Eliminate pass done" << endl;
            }

            // Reset member variables
            node->marked = false;
            node->relinearize = false;
            node->changedColStructure.clear();
            node->changedColContribution.clear();
            node->newColContribution.clear();
            node->reconstructCols.clear();
            node->editCols.clear();
            node->adoptedSkip = false;
            node->adoptedCols.clear();
        }
        
    }

}
*/

void CholeskyEliminationTree::print(std::ostream& os) {

    // map<size_t, Key> orderingToKey_;
    // for(int k = 0; k < ordering_.size(); k++) {
    //     orderingToKey_.insert({ordering_[k], k});
    // }
    // os << "Ordering: " << endl;
    // for(const auto p : orderingToKey_) {
    //     os << p.second << " ";
    // }
    // os << endl << endl;
    // os << "Jacobian:" << endl;
    // jacobian_.print(os);
    // os << "Jacobian end" << endl << endl;
    // os << "Cholesky:" << endl;
    // cholesky_.print(os);
    // os << "Cholesky end" << endl << endl;
    // // os << "b:" << endl;
    // // b_.print(os);
    // // os << "b end" << endl << endl;
    // // os << "y:" << endl;
    // // y_.print(os);
    // // os << "y end" << endl << endl;
    // os << "delta:" << endl;
    // delta_.print(os);
    // os << "delta end" << endl << endl;

    os << "Delta = Vector Values with " << orderingToKey_.size() << " elements" << endl;
    for(size_t i = 0; i < orderingToKey_.size(); i++) {
        os << i << ": " << delta_.block(i).transpose() << endl;
    }
    os << endl;
}

/*
void CholeskyEliminationTree::eliminateColumn(SparseColumnBlockMatrix* column_ptr, 
                                              Eigen::VectorXd* diagonalY) {
    SparseColumnBlockMatrix& column = *column_ptr;

    // cout << "new column" << endl;
    // column.print(cout);
    Block D = column.diagonalBlock();
    Eigen::LLT<Eigen::Ref<Block>, Eigen::Lower> llt(D);
    auto L = D.triangularView<Eigen::Lower>();

    L.solveInPlace(*diagonalY); // Ly = Atb

    if(column.height() > column.width()) {
        Block B = column.belowDiagonalBlocks();
        L.solveInPlace(B.transpose());
        // cout << "After solve in place" << endl;
        // column.print(cout);
        const vector<pair<Key, RowHeightPair>>& blockStartVec = column.blockStartVec();
        size_t remainingHeight = column.height() - column.width();

        // // 03/09/2023: Redo this part to do outer product at once
        vector<RowHeightPair> destBlockStart;
        RowMajorMatrix augmentedB(remainingHeight + 1, column.width());     // This is the augmented matrix with y
        augmentedB.block(0, 0, remainingHeight, column.width()) = B;
        augmentedB.block(remainingHeight, 0, 1, column.width()) = diagonalY->transpose();
        RowMajorMatrix outerProduct(remainingHeight + 1, remainingHeight + 1);
        outerProduct.setZero();
        outerProduct.selfadjointView<Eigen::Lower>().rankUpdate(augmentedB);

        for(size_t i = 1; i < blockStartVec.size(); i++) {   // start from block 1 instead of 0
            const Key otherKey = blockStartVec[i].first;
            auto& destColumn = cholesky_.column(otherKey);

            // Find the destination indices first
            destBlockStart.clear();
            size_t lastDestRow = INT_MAX;
            const auto& otherBlockStartMap = destColumn.blockStartMap();
            for(size_t j = i; j < blockStartVec.size(); j++) {
                const Key destKey = blockStartVec[j].first;
                const RowHeightPair p = otherBlockStartMap.at(destKey);
                const size_t destRow = p.first;
                const size_t destHeight = p.second;
                if(lastDestRow != destRow) {
                    destBlockStart.push_back({destRow, destHeight});
                }
                else {
                    // Merge two blocks together
                    destBlockStart.back().second += destHeight;
                }
                lastDestRow = destRow + destHeight;
            }

            // Take width off since we're not counting diagonal blocks
            const size_t otherRow = blockStartVec[i].second.first - column.width();
            const size_t otherHeight = blockStartVec[i].second.second;
            size_t curRow = otherRow;
            for(const RowHeightPair p : destBlockStart) {
                const size_t destRow = p.first;
                const size_t destHeight = p.second;
                destColumn.blockRange(destRow, destHeight) 
                    -= outerProduct.block(curRow, otherRow, destHeight, otherHeight);
                curRow += destHeight;
            }

        }
        
        // Update Atb, separate from the previous loop because the destinations could 
        // be out of order, and we don't want to update too many times
        const size_t lastRow = remainingHeight;
        size_t curCol = 0;
        for(size_t i = 1; i < blockStartVec.size(); i++) {
            const Key destKey = blockStartVec[i].first;
            const size_t destHeight = blockStartVec[i].second.second;
            y_.block(destKey) -= outerProduct.block(lastRow, curCol, 1, destHeight).transpose();
            curCol += destHeight;
        }

    }

    cout << "After eliminate column" << endl;
    column.print(cout);
}
*/

/*
void CholeskyEliminationTree::handleEdits(sharedNode node) {
    const Key key = node->key;
    SparseColumnBlockMatrix& column = cholesky_.column(key);
    // handle edits only if matrix is not new
    // If node is marked and there is a reordering
    // All ancestors will also be reconstructs, so don't have to worry about it
    if(column.matrixHeight() != 0) {
        assert(column.allocated());
        const auto& blockStartVec = column.blockStartVec();
        auto editIter = node->editCols.rbegin();  // editCols is in reverse order because we're going from the top to bot of the tree
        auto editIterEnd = node->editCols.rend();
        vector<RowHeightPair> destBlockStart(blockStartVec.size());
        // Scratchpad space for outer product
        RowMajorMatrix outerProductColumn(column.height(), column.width());

        for(size_t i = 1; i < blockStartVec.size(); i++) {
            const Key otherKey = blockStartVec[i].first;
            if(editIter != editIterEnd 
                    && *editIter == otherKey) {
                // Found ancestor key in reconstructCols
                editIter++;
                // cout << "Edit column " << otherKey << endl;
            }
            else {
                continue;
            }
            const size_t otherRow = blockStartVec[i].second.first;
            const size_t otherHeight = blockStartVec[i].second.second;

            auto& destColumn = cholesky_.column(otherKey);

            // cout << "In reconstruct. orig est column = " << endl;
            // destColumn.print(cout);

            // Find the destination indices first
            const auto& otherBlockStartMap = destColumn.blockStartMap();
            for(size_t j = i; j < blockStartVec.size(); j++) {
                const Key destKey = blockStartVec[j].first;
                destBlockStart[j] = otherBlockStartMap.at(destKey);
                assert(column.blockStartMap().at(destKey).second == destBlockStart[j].second);
            }

            size_t remainingHeight = column.height() - otherRow;
            Block curBlock = column.blockRange(otherRow, otherHeight);
            Block remainingBlocks = column.blockRange(otherRow, remainingHeight);

            outerProductColumn.block(otherRow, 0, 
                    remainingHeight, column.width()).noalias() 
                = remainingBlocks * curBlock.transpose();


            size_t curRow = otherRow;
            for(size_t j = i; j < blockStartVec.size(); j++) {
                const size_t destRow = destBlockStart[j].first;
                const size_t destHeight = destBlockStart[j].second;
                // For edit, we add the block back
                destColumn.blockRange(destRow, destHeight) 
                    += outerProductColumn.block(curRow, 0, destHeight, otherHeight);
                curRow += destHeight;

            }
        }

        // FIXME: handle edit needs to fix Atb
        if(column.hasBelowDiagonalBlocks() && !node->editCols.empty()) {
            editIter = node->editCols.rbegin();  
            editIterEnd = node->editCols.rend();

            // Edit needs to fix Atb. check if in editCols
            Eigen::VectorXd outerprodRhs = column.belowDiagonalBlocks() * y_.block(key);

            size_t curRow = 0;
            for(size_t i = 1; i < blockStartVec.size(); i++) {
                const size_t destKey = blockStartVec[i].first;
                const size_t destHeight = blockStartVec[i].second.second;

                if(editIter != editIterEnd
                        && *editIter == destKey) {
                    // Found ancestor key in editIter
                    editIter++;
                    // ADD outerProduct
                    y_.block(destKey) += outerprodRhs.block(curRow, 0, destHeight, 1);
                }
                curRow += destHeight;
            }
        }
    }
}

void CholeskyEliminationTree::handleReconstructs(sharedNode node) {
    const Key key = node->key;
    auto& column = cholesky_.column(key);
    const auto& blockStartVec = column.blockStartVec();
    auto reconstructIter = node->reconstructCols.rbegin();  // reconstructCols is in reverse order because we're going from the top to bot of the tree
    auto reconstructIterEnd = node->reconstructCols.rend();
    vector<RowHeightPair> destBlockStart(blockStartVec.size());
    // Scratchpad space for outer product
    RowMajorMatrix outerProductColumn(column.height(), column.width());
    for(size_t i = 1; i < blockStartVec.size(); i++) {
        const Key otherKey = blockStartVec[i].first;
        if(reconstructIter != reconstructIterEnd 
                && *reconstructIter == otherKey) {
            // Found ancestor key in reconstructCols
            reconstructIter++;
            // cout << "Reconstruct column " << otherKey << endl;
        }
        else {
            continue;
        }
        const size_t otherRow = blockStartVec[i].second.first;
        const size_t otherHeight = blockStartVec[i].second.second;

        auto& destColumn = cholesky_.column(otherKey);

        // Find the destination indices first
        const auto& destBlockStartMap = destColumn.blockStartMap();
        for(size_t j = i; j < blockStartVec.size(); j++) {
            const Key destKey = blockStartVec[j].first;
            destBlockStart[j] = otherBlockStartMap.at(destKey);
            assert(column.blockStartMap().at(destKey).second == destBlockStart[j].second);
        }

        size_t remainingHeight = column.height() - otherRow;
        Block curBlock = column.blockRange(otherRow, otherHeight);
        Block remainingBlocks = column.blockRange(otherRow, remainingHeight);

        outerProductColumn.block(otherRow, 0, 
                remainingHeight, column.width()).noalias()
            = remainingBlocks * curBlock.transpose();


        size_t curRow = otherRow;
        for(size_t j = i; j < blockStartVec.size(); j++) {
            const size_t destRow = destBlockStart[j].first;
            const size_t destHeight = destBlockStart[j].second;
            destColumn.blockRange(destRow, destHeight) 
                -= outerProductColumn.block(curRow, 0, destHeight, otherHeight);
            curRow += destHeight;

        }
    }
    // TODO: Add Atb as a block in the column
    if(column.hasBelowDiagonalBlocks() && !node->reconstructCols.empty()) {
        reconstructIter = node->reconstructCols.rbegin();  
        reconstructIterEnd = node->reconstructCols.rend();

        // Reconstruct needs to fix Atb. FIXME: check if in reconstruct Cols
        Eigen::VectorXd outerprodRhs = column.belowDiagonalBlocks() * y_.block(key);

        size_t curRow = 0;
        for(size_t i = 1; i < blockStartVec.size(); i++) {
            const size_t destKey = blockStartVec[i].first;
            const size_t destHeight = blockStartVec[i].second.second;

            if(reconstructIter != reconstructIterEnd 
                    && *reconstructIter == destKey) {
                // Found ancestor key in reconstructCols
                reconstructIter++;
                y_.block(destKey) -= outerprodRhs.block(curRow, 0, destHeight, 1);
            }
            curRow += destHeight;
        }
    }
}
*/

bool CholeskyEliminationTree::gatherColumns(
        sharedClique clique, 
        ColMajorMatrix* m,
        size_t* totalWidth,
        size_t* totalHeight,
        vector<Key>* keys,
        vector<size_t>* widths,
        vector<size_t>* heights,
        BlockIndexVector* blockIndices,
        vector<LowerTriangularColumnMatrix*>* columns) {

    // cout << "gatherColumns" << endl;

    // TODO: Optimize this for when clique only has one column
    // At this point, some new columns may not be allocated nor fully allocated
    *totalWidth = 0;
    for(sharedNode node : clique->nodes) {
        auto column_ptr = &cholesky_.column(node->key);
        if(!column_ptr->allocated()) {
            // cout << "skipped " << node->key << endl;
            assert(node == clique->back());
            break;
        }
        else {
            // cout << "gathering " << node->key << endl;
        }
        keys->push_back(node->key);
        columns->push_back(column_ptr);
        widths->push_back(column_ptr->width());
        heights->push_back(column_ptr->matrixHeight());
        *totalWidth += widths->back();
    }

    if(columns->empty()) {
        // Generally the final column is in a clique with the second
        // to last column
        return false;
    }

    *totalHeight = heights->front();   // This needs to account for the b row

    *m = ColMajorMatrix(*totalHeight, *totalWidth);

    size_t curCol = 0, curRow = 0;
    for(size_t i = 0; i < keys->size(); i++) {
        m->block(curRow, curCol, (*heights)[i], (*widths)[i]) 
            = (*columns)[i]->blockRange(0, (*heights)[i]);

        curCol += (*widths)[i];
        curRow += (*widths)[i];
    }

    // Need to zero out the upper diagonal
    Block Dblock = m->block(0, 0, *totalWidth, *totalWidth);
    Eigen::VectorXd D = Dblock.diagonal();
    Dblock.triangularView<Eigen::Upper>().setZero();
    Dblock.diagonal() = D;

    // Populate blockIndices for the new column, it can just be diagonal block
    // plus the rest of the last column's blockIndices
    *blockIndices = vector<pair<Key, RowHeightPair>>(1, {keys->front(), {0, *totalWidth}});
    const auto& lastBlockIndices = columns->back()->blockIndices();
    curRow = *totalWidth;
    for(size_t i = 1; i < lastBlockIndices.size(); i++) {
        const auto& p = lastBlockIndices[i];
        const Key key = p.first;
        const size_t height = p.second.second;
        blockIndices->push_back({key, {curRow, height}});
        curRow += height;
    }

    return true;
}

void CholeskyEliminationTree::scatterColumns(const ColMajorMatrix& m, 
        const std::vector<size_t>& widths,
        const std::vector<size_t>& heights,
        const std::vector<LowerTriangularColumnMatrix*>& columns) {
    // cout << "In scatterColumns" << endl;
    size_t curRow = 0;
    for(size_t i = 0; i < widths.size(); i++) {
        auto& destColumn = *(columns[i]);
        size_t width = widths[i];
        size_t height = heights[i];
        // cout << "key = " << destColumn.key() << " width = " << width << " height = " << height << " col height = " << destColumn.height() << endl;
        destColumn.blockRange(0, height) 
            = m.block(curRow, curRow, height, width);
        curRow += width;
    }
}

void CholeskyEliminationTree::editAndRestoreFromClique(sharedClique clique) {
    // cout << "in edit and restore. editCols.size = " << clique->editCols.size() << " " << clique->is_reconstruct << endl;
    vector<Key> editCols;
    editCols.reserve(clique->back()->colStructure.size());
    for(Key key : clique->back()->colStructure) {
        if(markedStatus_[key] == EDIT) {
            editCols.push_back(key);
        }
    }

    // cout << "editCols: ";
    // for(Key k : editCols) {
    //     cout << k << " ";
    // }
    // cout << endl;

    if(editCols.empty() && clique->is_reconstruct) {
        // If no edit cols and no need to restore column
        return;
    }

    // First gather all the columns
    const Key firstKey = clique->front()->key;
    size_t totalWidth = 0;
    size_t totalHeight = 0;
    vector<Key> keys;
    vector<size_t> widths;
    vector<size_t> heights;
    BlockIndexVector blockStartVec;
    vector<LowerTriangularColumnMatrix*> columns;

    ColMajorMatrix m;

    // Here we do need to care if cannot gather
    bool gathered = gatherColumns(clique, &m, &totalWidth, &totalHeight, 
                                  &keys, &widths, &heights, &blockStartVec, &columns);
    // cout << "edit and restore 1" << endl;

    if(!gathered) {
        // cout << "New clique. nothing to Gather" << endl;
        // exit(1);
        return;
    }

    if(!editCols.empty()) {
        editOrReconstructFromColumn(m, blockStartVec, editCols, keys, 1);
        // editFromColumn(m, blockStartVec, editCols, keys);
    }
    if(!clique->is_reconstruct) {
        restoreColumn(m, totalWidth, totalHeight);
        scatterColumns(m, widths, heights, columns);
    }
}

void CholeskyEliminationTree::restoreColumn(ColMajorMatrix& m, 
        size_t totalWidth, size_t totalHeight) {

    Block D = m.block(0, 0, totalWidth, totalWidth);
    Block col = m.block(0, 0, totalHeight, totalWidth);
    // Reset D's upper triangular matrix
    Eigen::VectorXd d = D.diagonal();
    D.triangularView<Eigen::Upper>().setZero();
    D.diagonal() = d;

    // ColMajorMatrix colCopy = col;
    // ColMajorMatrix DTCopy = D.transpose();

    // We did L^-1 B^T before, now we want to do L*B^T which is B*L^T
    // col.noalias() = colCopy * DTCopy;
    col *= D.transpose();

}

void CholeskyEliminationTree::editOrReconstructFromColumn(
    const ColMajorMatrix& m,
    const vector<pair<Key, RowHeightPair>>& blockStartVec,
    const vector<Key>& destCols,
    const vector<Key>& gatheredKeys,
    double sign) {

    bool processGrouped = true;
    size_t destSize = destCols.size();
    size_t totalSize = blockStartVec.size() - 1;
    if(destSize * 2 < totalSize) {
        processGrouped = false;
        // cout << "edit single" << endl;
    }

    // Maybe we should have 2 ways of doing this, one where there are not
    // a lot of edit keys, one where there are
    size_t firstRow = blockStartVec[1].second.first;
    size_t bHeight = m.rows() - firstRow;
    size_t bWidth = blockStartVec[0].second.second;             // diagonal height

    auto b = m.block(firstRow, 0, bHeight, bWidth);
    ColMajorMatrix outerProduct;
    if(processGrouped) {
        outerProduct = ColMajorMatrix(bHeight, bHeight);
        outerProduct.setZero();
        outerProduct.selfadjointView<Eigen::Lower>().rankUpdate(b, sign);
    }

    auto destIter = destCols.begin();
    auto destEnd = destCols.end();
    for(size_t i = 1; i < blockStartVec.size(); i++) {
        // Start from index 1 to skip diagonal block
        const Key otherKey = blockStartVec[i].first;
        if(destIter != destEnd
                && *destIter == otherKey) {
            // Found ancestor key 
            destIter++;
        }
        else {
            continue;
        }

        const auto p = blockStartVec[i].second;
        const size_t otherRow = p.first;
        const size_t otherHeight = p.second;

        auto& destColumn = cholesky_.column(otherKey);
    
        // Find destination indices first
        vector<RowHeightPair> destBlockStart;
        const auto& destBlockIndices = destColumn.blockIndices();
        findDestBlocks(blockStartVec.begin() + i,
                       blockStartVec.end(),
                       destBlockIndices,
                       &destBlockStart);

        // Start from row 0 for edit since we don't edit the diagonal
        size_t curRow = 0;
        size_t curCol = otherRow - bWidth;

        if(!processGrouped) {
            outerProduct = sign * m.block(otherRow, 0, bHeight - curCol, bWidth) * m.block(otherRow, 0, otherHeight, bWidth).transpose();
        }

        Block outerProductColumn = processGrouped? outerProduct.block(curCol, curCol, bHeight - curCol, otherHeight) : outerProduct.block(0, 0, bHeight - curCol, otherHeight);
        // if(editGrouped) {
        //     outerProductColumn 
        //         = outerProduct.block(curCol, curCol, bHeight - curCol, otherHeight);
        // }

        for(const auto& p : destBlockStart) {
            const size_t destRow = p.first;
            const size_t destHeight = p.second;
            // For edit, we ADD the block 
            // destColumn.blockRange(destRow, destHeight)
            //     += outerProduct.block(curRow, curCol, destHeight, otherHeight);
            destColumn.blockRange(destRow, destHeight)
                += outerProductColumn.block(curRow, 0, destHeight, otherHeight);
            curRow += destHeight;
        }
    }
}

void CholeskyEliminationTree::editFromColumn(const ColMajorMatrix& m,
    const vector<pair<Key, RowHeightPair>>& blockStartVec,
    const vector<Key>& editCols,
    const vector<Key>& gatheredKeys) {
    assert(!editCols.empty());
    // assert(editCols.back() != LAST_ROW);

    assert(blockStartVec.size() >= 2);
    assert(sorted_no_duplicates(editCols));

    bool editGrouped = true;
    size_t editSize = editCols.size();
    size_t totalSize = blockStartVec.size() - 1;
    if(editSize * 6 < totalSize) {
        editGrouped = false;
        // cout << "edit single" << endl;
    }

    // Maybe we should have 2 ways of doing this, one where there are not
    // a lot of edit keys, one where there are
    size_t firstRow = blockStartVec[1].second.first;
    size_t bHeight = m.rows() - firstRow;
    size_t bWidth = blockStartVec[0].second.second;             // diagonal height

    auto b = m.block(firstRow, 0, bHeight, bWidth);
    ColMajorMatrix outerProduct;
    if(editGrouped) {
        outerProduct = ColMajorMatrix(bHeight, bHeight);
        outerProduct.setZero();
        outerProduct.selfadjointView<Eigen::Lower>().rankUpdate(b);
    }

    auto editIter = editCols.begin();
    auto editEnd = editCols.end();
    for(size_t i = 1; i < blockStartVec.size(); i++) {
        // Start from index 1 to skip diagonal block
        const Key otherKey = blockStartVec[i].first;
        if(editIter != editEnd
                && *editIter == otherKey) {
            // Found ancestor key 
            editIter++;
        }
        else {
            continue;
        }

        const auto p = blockStartVec[i].second;
        const size_t otherRow = p.first;
        const size_t otherHeight = p.second;

        auto& destColumn = cholesky_.column(otherKey);
    
        // Find destination indices first
        vector<RowHeightPair> destBlockStart;
        const auto& destBlockIndices = destColumn.blockIndices();
        findDestBlocks(blockStartVec.begin() + i,
                       blockStartVec.end(),
                       destBlockIndices,
                       &destBlockStart);

        // Start from row 0 for edit since we don't edit the diagonal
        size_t curRow = 0;
        size_t curCol = otherRow - bWidth;

        if(!editGrouped) {
            outerProduct = m.block(otherRow, 0, bHeight - curCol, bWidth) * m.block(otherRow, 0, otherHeight, bWidth).transpose();
        }

        Block outerProductColumn = editGrouped? outerProduct.block(curCol, curCol, bHeight - curCol, otherHeight) : outerProduct.block(0, 0, bHeight - curCol, otherHeight);
        // if(editGrouped) {
        //     outerProductColumn 
        //         = outerProduct.block(curCol, curCol, bHeight - curCol, otherHeight);
        // }

        for(const auto& p : destBlockStart) {
            const size_t destRow = p.first;
            const size_t destHeight = p.second;
            // For edit, we ADD the block 
            // destColumn.blockRange(destRow, destHeight)
            //     += outerProduct.block(curRow, curCol, destHeight, otherHeight);
            destColumn.blockRange(destRow, destHeight)
                += outerProductColumn.block(curRow, 0, destHeight, otherHeight);
            curRow += destHeight;
        }

        // // DEBUG
        // for(const Key k : gatheredKeys) {
        //     assert(descendantStatus[otherKey].find(k) != descendantStatus[otherKey].end());
        //     assert(descendantStatus[otherKey].at(k) == EDIT);
        //     descendantStatus[otherKey].at(k) = READY;
        // }
    }
}

/*
void CholeskyEliminationTree::editFromColumn(const ColMajorMatrix& m,
    const vector<pair<Key, RowHeightPair>>& blockStartVec,
    const vector<Key>& editCols,
    const vector<Key>& gatheredKeys) {
    assert(!editCols.empty());
    // assert(editCols.back() != LAST_ROW);

    assert(blockStartVec.size() >= 2);
    assert(sorted_no_duplicates(editCols));

    bool editGrouped = true;
    size_t editSize = editCols.size();
    size_t totalSize = blockStartVec.size() - 1;
    if(editSize * 2 < totalSize) {
        editGrouped = false;
    }

    // Maybe we should have 2 ways of doing this, one where there are not
    // a lot of edit keys, one where there are
    size_t firstRow = blockStartVec[1].second.first;
    size_t bHeight = m.rows() - firstRow;
    size_t bWidth = blockStartVec[0].second.second;             // diagonal height

    auto b = m.block(firstRow, 0, bHeight, bWidth);
    ColMajorMatrix outerProduct(bHeight, bHeight);
    outerProduct.setZero();

    outerProduct.selfadjointView<Eigen::Lower>().rankUpdate(b);

    auto editIter = editCols.begin();
    auto editEnd = editCols.end();
    for(size_t i = 1; i < blockStartVec.size(); i++) {
        // Start from index 1 to skip diagonal block
        const Key otherKey = blockStartVec[i].first;
        if(editIter != editEnd
                && *editIter == otherKey) {
            // Found ancestor key 
            editIter++;
        }
        else {
            continue;
        }

        const auto p = blockStartVec[i].second;
        const size_t otherRow = p.first;
        const size_t otherHeight = p.second;

        auto& destColumn = cholesky_.column(otherKey);
    
        // Find destination indices first
        vector<RowHeightPair> destBlockStart;
        const auto& destBlockIndices = destColumn.blockIndices();
        findDestBlocks(blockStartVec.begin() + i,
                       blockStartVec.end(),
                       destBlockIndices,
                       &destBlockStart);

        // Start from row 0 for edit since we don't edit the diagonal
        size_t curRow = otherRow - bWidth;
        size_t curCol = otherRow - bWidth;
        for(const auto& p : destBlockStart) {
            const size_t destRow = p.first;
            const size_t destHeight = p.second;
            // For edit, we ADD the block 
            destColumn.blockRange(destRow, destHeight)
                += outerProduct.block(curRow, curCol, destHeight, otherHeight);
            curRow += destHeight;
        }

        // // DEBUG
        // for(const Key k : gatheredKeys) {
        //     assert(descendantStatus[otherKey].find(k) != descendantStatus[otherKey].end());
        //     assert(descendantStatus[otherKey].at(k) == EDIT);
        //     descendantStatus[otherKey].at(k) = READY;
        // }
    }
}
*/

bool CholeskyEliminationTree::reconstructFromClique(sharedClique clique) {
    // First gather all the columns
    // then call reconstruct column
    
    const Key firstKey = clique->front()->key;
    size_t totalWidth = 0;
    size_t totalHeight = 0;
    vector<Key> keys;
    vector<size_t> widths;
    vector<size_t> heights;
    BlockIndexVector blockStartVec;
    vector<LowerTriangularColumnMatrix*> columns;

    vector<Key> reconstructCols;
    reconstructCols.reserve(clique->back()->colStructure.size());
    for(Key key : clique->back()->colStructure) {
        if(markedStatus_[key] == RECONSTRUCT) {
            // cout << "Pushed " << key << " to reconstructCols" << endl;
            reconstructCols.push_back(key);
        }
    }
    assert(sorted_no_duplicates(clique->back()->colStructure));

    // cout << "reconstructCols: ";
    // for(Key k : reconstructCols) {
    //     cout << "[" << k << " " << ordering_[k] << "] ";
    // }
    // cout << endl;

    assert(sorted_no_duplicates(reconstructCols));

    if(reconstructCols.empty()) {
        return false;
    }

    ColMajorMatrix m;

    bool gathered = gatherColumns(clique, &m, &totalWidth, &totalHeight, 
                                  &keys, &widths, &heights, &blockStartVec, &columns);
    assert(gathered);

    editOrReconstructFromColumn(m, blockStartVec, reconstructCols, keys, -1);
    // reconstructFromColumn(m, blockStartVec, reconstructCols, keys);

    return true;
}

void CholeskyEliminationTree::reconstructFromColumn(const ColMajorMatrix& m,
    const vector<std::pair<Key, RowHeightPair>>& blockStartVec,
    const vector<Key>& reconstructCols,
    const vector<Key>& gatheredKeys) {
    assert(!reconstructCols.empty());
    assert(sorted_no_duplicates(reconstructCols));
    // assert(reconstructCols.back() != LAST_ROW);

    assert(blockStartVec.size() >= 2);

    bool reconstructGrouped = true;
    size_t reconstructSize = reconstructCols.size();
    size_t totalSize = blockStartVec.size() - 1;
    if(reconstructSize * 6 < totalSize) {
        reconstructGrouped = false;
        // cout << "edit single" << endl;
    }

    // Maybe we should have 2 ways of doing this, one where there are not
    // a lot of reconstruct keys, one where there are
    size_t firstRow = blockStartVec[1].second.first;
    size_t bHeight = m.rows() - firstRow;
    size_t bWidth = blockStartVec[0].second.second;             // diagonal height

    auto b = m.block(firstRow, 0, bHeight, bWidth);
    ColMajorMatrix outerProduct;
    if(reconstructGrouped) {
        outerProduct = ColMajorMatrix(bHeight, bHeight);
        outerProduct.setZero();
        outerProduct.selfadjointView<Eigen::Lower>().rankUpdate(b);
    }

    auto reconstructIter = reconstructCols.begin();
    auto reconstructEnd = reconstructCols.end();
    for(size_t i = 1; i < blockStartVec.size(); i++) {
        // Start from index 1 to skip diagonal block
        const Key otherKey = blockStartVec[i].first;
        if(reconstructIter != reconstructEnd
                && *reconstructIter == otherKey) {
            // Found ancestor key 
            reconstructIter++;
            // cout << "Reconstruct column " << otherKey << endl;
        }
        else {
            continue;
        }

        const auto p = blockStartVec[i].second;
        const size_t otherRow = p.first;
        const size_t otherHeight = p.second;

        auto& destColumn = cholesky_.column(otherKey);
    
        // Find destination indices first
        vector<RowHeightPair> destBlockStart;
        const auto& destBlockIndices = destColumn.blockIndices();
        findDestBlocks(blockStartVec.begin() + i,
                       blockStartVec.end(),
                       destBlockIndices,
                       &destBlockStart);

        // Start from row 0 for reconstruct since we don't reconstruct the diagonal
        size_t curRow = 0;
        size_t curCol = otherRow - bWidth;

        if(!reconstructGrouped) {
            outerProduct = m.block(otherRow, 0, bHeight - curCol, bWidth) * m.block(otherRow, 0, otherHeight, bWidth).transpose();
        }

        Block outerProductColumn = reconstructGrouped? outerProduct.block(curCol, curCol, bHeight - curCol, otherHeight) : outerProduct.block(0, 0, bHeight - curCol, otherHeight);

        // size_t curRow = otherRow - bWidth;
        // size_t curCol = otherRow - bWidth;
        for(const auto& p : destBlockStart) {
            const size_t destRow = p.first;
            const size_t destHeight = p.second;
            // cout << "destKey = " << destColumn.key() << " destRow = " << destRow << " destHeight = " << destHeight << endl;
            // For reconstruct, we SUBTRACT the block 
            // destColumn.blockRange(destRow, destHeight)
            //     -= outerProduct.block(curRow, curCol, destHeight, otherHeight);
            destColumn.blockRange(destRow, destHeight)
                -= outerProductColumn.block(curRow, 0, destHeight, otherHeight);
            curRow += destHeight;
        }
    }
}

void CholeskyEliminationTree::eliminateClique(sharedClique clique) {
    // cout << "CholeskyEliminationTree::eliminateClique()" << endl;
    // First gather all the columns
    // then call eliminate column
    const Key firstKey = clique->front()->key;
    size_t totalWidth = 0;
    size_t totalHeight = 0;
    vector<Key> keys;
    vector<size_t> widths;
    vector<size_t> heights;
    BlockIndexVector blockStartVec;
    vector<LowerTriangularColumnMatrix*> columns;

    ColMajorMatrix m;

    gatherColumns(clique, &m, &totalWidth, &totalHeight, 
                  &keys, &widths, &heights, &blockStartVec, &columns);

    eliminateColumn(m, totalWidth, totalHeight, blockStartVec, keys);

    scatterColumns(m, widths, heights, columns);
}

void CholeskyEliminationTree::eliminateColumn(ColMajorMatrix& m,
        const size_t totalWidth,
        const size_t totalHeight,
        const vector<std::pair<Key, RowHeightPair>>& blockStartVec,
        const vector<Key>& gatheredKeys) {

    // cout << "Before eliminate column. m = " << endl << m << endl << endl;  

    size_t diagonalWidth = totalWidth;
    Block D = m.block(0, 0, diagonalWidth, diagonalWidth);
    Eigen::LLT<Eigen::Ref<Block>, Eigen::Lower> llt(D);
    if(llt.info() == Eigen::NumericalIssue) {
        cout << "Diagonal block not positive definite!" << endl;
        Key firstKey = blockStartVec[0].first;
        cout << "First key: " << firstKey << " clique is_reconstruct? " << nodes_[firstKey]->clique->is_reconstruct << endl;
        exit(1);
    }
    auto L = D.triangularView<Eigen::Lower>();

    if(totalHeight > totalWidth) {  
        size_t remainingHeight = totalHeight - totalWidth;
        Block B = m.block(totalWidth, 0, remainingHeight, totalWidth);
        L.solveInPlace(B.transpose());

        // if(nodes_.size() >= 2500) {
        //     cout << remainingHeight << " " << totalWidth << endl;
        // }

        ColMajorMatrix outerProduct(remainingHeight, remainingHeight);
        outerProduct.setZero();
        outerProduct.selfadjointView<Eigen::Lower>().rankUpdate(B);

        for(size_t i = 1; i < blockStartVec.size() - 1; i++) {
            // Start from index 1 to skip diagonal block
            // end at -1 as the b row cannot be a source of outerproduct
            const Key otherKey = blockStartVec[i].first;

            const auto p = blockStartVec[i].second;
            const size_t otherRow = p.first;
            const size_t otherHeight = p.second;

            auto& destColumn = cholesky_.column(otherKey);

            // Find destination indices first
            vector<RowHeightPair> destBlockStart;
            const auto& destBlockIndices = destColumn.blockIndices();
            // cout << "eliminate column findDestBlocks otherKey = " << otherKey << endl; 
            findDestBlocks(blockStartVec.begin() + i,
                    blockStartVec.end(),
                    destBlockIndices,
                    &destBlockStart);

            // Start from row 0 for reconstruct since we don't reconstruct the diagonal
            size_t curRow = otherRow - diagonalWidth;
            size_t curCol = otherRow - diagonalWidth;
            for(const auto& p : destBlockStart) {
                const size_t destRow = p.first;
                const size_t destHeight = p.second;
                // cout << "destColumn " << destColumn.key() << " destRow = " << destRow 
                //      << " destHeight = " << destHeight << endl << endl;
                // For eliminate, we SUBTRACT the block 
                destColumn.blockRange(destRow, destHeight)
                    -= outerProduct.block(curRow, curCol, destHeight, otherHeight);
                curRow += destHeight;
            }
        }
    }
    // cout << "After eliminate column. m = " << endl << m << endl << endl;  

}

template<typename Iterator>
void CholeskyEliminationTree::findDestBlocks(Iterator start, Iterator end, 
        const BlockIndexVector& destBlockIndices,
        std::vector<RowHeightPair>* destBlockStart) {

    destBlockStart->clear();

    auto destIt = destBlockIndices.begin();
    auto destEnd = destBlockIndices.end();
    size_t lastRow = LAST_ROW - 1;    // Use to merge destination blocks
    for(; start != end; start++) {
        const Key destKey = start->first;

        while(destIt->first != destKey) {
            destIt++;
            assert(destIt != destEnd);
        }

        size_t destRow = destIt->second.first;
        size_t destHeight = destIt->second.second;

        if(destRow != lastRow) {
            destBlockStart->push_back({destRow, destHeight});
        }
        else {
            // merge blocks
            destBlockStart->back().second += destHeight;
        }
        lastRow = destRow + destHeight;
    }
}

void CholeskyEliminationTree::simpleAtB(Block A, Block B, Block C, 
                                        size_t m, size_t k, size_t n, bool add) {
    if(add) {
        for(size_t i = 0; i < m; i++) {
            for(size_t j = 0; j < n; j++) {
                for(size_t t = 0; t < k; t++) {
                    C(i, j) += A(t, i) * B(t, j);
                }
            }
        }
    }
}

// subtract AtA from each column of the clique
void CholeskyEliminationTree::prepareEditClique(sharedClique clique) {
    for(sharedNode node : clique->nodes) {
        prepareEditColumn(node);
    }
}

void CholeskyEliminationTree::prepareEditColumn(sharedNode node) {
    const Key key = node->key;
    // cout << "[CholeskyEliminationTree] prepareEditColumn. key = " << key << endl;
    LowerTriangularColumnMatrix& column = cholesky_.column(key);

    // Sort rows in order
    map<Key, vector<FactorIndex>, OrderingLess> destIndices(orderingLess);

    for(const FactorIndex factorIndex : node->factorIndices) {
        // Unlinearized factors don't need to be reset
        if(factorLinearizeStatus_[factorIndex] == RELINEARIZED) {
            sharedFactor factor = factors_[factorIndex];
            for(const Key otherKey : factor->keys()) {
                // This key need to be added to AtA
                if(!orderingLess(otherKey, key)) {
                    destIndices[otherKey].push_back(factorIndex);
                }
            }

            // Do Atb here since it only needs to be done once per factor
            // SUBTRACT Atb for edits
            column.lastRow()
                -= jacobian_.block(factorIndex, -1).transpose() * jacobian_.block(factorIndex, key);
        }
    }

    auto it = column.blockIndices().begin();
    auto itEnd = column.blockIndices().end();

    // Now subtract factors
    for(auto& p : destIndices) {
        Key otherKey = p.first;
        while(it->first != otherKey) {
            it++;
            assert(it != itEnd);    // Our block index must already have allocated the block
        }
        const size_t& destRow = it->second.first;
        const size_t& destHeight = it->second.second;

        for(FactorIndex factorIndex : p.second) {
            // SUBTRACT AtA block
            column.blockRange(destRow, destHeight)
                -= jacobian_.block(factorIndex, otherKey).transpose() * jacobian_.block(factorIndex, key);
        }
    }
}

void CholeskyEliminationTree::prepareEliminateClique(sharedClique clique, const Values& theta) {
    for(sharedNode node : clique->nodes) {
        prepareEliminateColumn(node, clique->is_reconstruct, theta);
    }
}

/*
void CholeskyEliminationTree::prepareEliminateColumn(sharedNode node, 
        bool is_reconstruct, const Values& theta) {
    // All factors of this node should be relinearized
    const Key key = node->key;
    LowerTriangularColumnMatrix& column = cholesky_.column(key);

    unordered_map<Key, RowHeightPair> blockIndexMap;
    for(auto& p : column.blockIndices()) {
        blockIndexMap.insert(p);
    }

    for(const FactorIndex factorIndex : node->factorIndices) {
        if(!is_reconstruct && factorLinearizeStatus_[factorIndex] == LINEARIZED) {
            // Can skip factor if node is edit and factor has already been linearized
            // But cannot skip if the factor is newly linearized, in that
            // case we have to add AtA
            continue;
        }

        sharedFactor factor = factors_[factorIndex];
        if(factorLinearizeStatus_[factorIndex] != LINEARIZED 
                && factorLinearizeStatus_[factorIndex] != NEWLINEARIZED) {
            vector<Matrix> A;
            Vector b(factor->dim());
            factor->linearizeToMatrix(theta, &A, &b);
            for(int i = 0; i < factor->size(); i++) {
                const Key factorKey = factor->keys()[i];
                jacobian_.block(factorIndex, factorKey) = A[i];
            }
            jacobian_.block(factorIndex, -1) = b;
            factorLinearizeStatus_[factorIndex] = NEWLINEARIZED;
        }
        bool lastKeyFlag = true;
        for(Key otherKey : factor->keys()) {
            if(!orderingLess(otherKey, key)) {
                auto& p = blockIndexMap.at(otherKey);
                size_t& destRow = p.first;
                size_t& destHeight = p.second;
                column.blockRange(destRow, destHeight)
                    += jacobian_.block(factorIndex, otherKey).transpose() * jacobian_.block(factorIndex, key);

                if(key != otherKey) {
                    lastKeyFlag = false;
                }
            }
        }
        // Atb row
        column.lastRow()
            += jacobian_.block(factorIndex, -1).transpose() * jacobian_.block(factorIndex, key);
        if(lastKeyFlag) {
            factorLinearizeStatus_[factorIndex] = LINEARIZED;
        }
    }
}
*/

void CholeskyEliminationTree::prepareEliminateColumn(sharedNode node, 
        bool is_reconstruct, const Values& theta) {
    // All factors of this node should be relinearized
    const Key key = node->key;
    LowerTriangularColumnMatrix& column = cholesky_.column(key);

    // Sort rows in order
    map<Key, vector<FactorIndex>, OrderingLess> destIndices(orderingLess);

    for(const FactorIndex factorIndex : node->factorIndices) {
        if(!is_reconstruct && factorLinearizeStatus_[factorIndex] == LINEARIZED) {
            // Can skip factor if node is edit and factor has already been linearized
            // But cannot skip if the factor is newly linearized, in that
            // case we have to add AtA
            continue;
        }
        sharedFactor factor = factors_[factorIndex];

        if(factorLinearizeStatus_[factorIndex] != LINEARIZED 
                && factorLinearizeStatus_[factorIndex] != NEWLINEARIZED) {
            vector<Matrix> A;
            Vector b(factor->dim());
            factor->linearizeToMatrix(theta, &A, &b);
            for(int i = 0; i < factor->size(); i++) {
                const Key factorKey = factor->keys()[i];
                jacobian_.block(factorIndex, factorKey) = A[i];
            }
            jacobian_.block(factorIndex, -1) = b;
            factorLinearizeStatus_[factorIndex] = NEWLINEARIZED;
        }

        bool lastKeyFlag = true;    // If this key is the last key in the factor
        for(const Key otherKey : factor->keys()) {
            // This key need to be added to AtA
            if(!orderingLess(otherKey, key)) {
                destIndices[otherKey].push_back(factorIndex);
                if(otherKey != key) {
                    // If there is a key higher order than this key
                    lastKeyFlag = false;
                }
            }
        }
        if(lastKeyFlag) {
            factorLinearizeStatus_[factorIndex] = LINEARIZED;
        }

        // Do Atb here since it only needs to be done once per factor
        // ADD Atb
        column.lastRow()
            += jacobian_.block(factorIndex, -1).transpose() * jacobian_.block(factorIndex, key);
        cout << "Atb" << key << endl << jacobian_.block(factorIndex, -1).transpose() * jacobian_.block(factorIndex, key) << endl << endl;
        cout << jacobian_.block(factorIndex, -1) << endl << endl << jacobian_.block(factorIndex, key) << endl << endl;
    }

    auto it = column.blockIndices().begin();
    auto itEnd = column.blockIndices().end();

    // Now ADD AtA blocks
    for(auto& p : destIndices) {
        Key otherKey = p.first;
        while(it->first != otherKey) {
            it++;
            assert(it != itEnd);    // Our block index must already have allocated the block
        }
        const size_t& destRow = it->second.first;
        const size_t& destHeight = it->second.second;

        for(FactorIndex factorIndex : p.second) {
            // ADD AtA block
            column.blockRange(destRow, destHeight) += jacobian_.block(factorIndex, otherKey).transpose() * jacobian_.block(factorIndex, key);
        }
    }
}

void CholeskyEliminationTree::reorderClique(sharedClique clique) {
    // cout << "in reorder clique" << endl;
    // This may not be true as we can have cliques that are newly marked
    // not be part of the reordering
    // assert(!clique->marked);

    for(sharedNode node : clique->nodes) {
        assert(node->ordering_version != ordering_version_);
        node->ordering_version = ordering_version_;
    }

    // Find the keys that have been reordered and the lowest reordered key
    auto& colStructure = clique->front()->colStructure; 
    vector<Key> oldColStructure = clique->front()->colStructure; // make a copy to compare

    std::sort(colStructure.begin(), colStructure.end(), orderingLess);

    for(size_t i = 1; i < clique->nodes.size(); i++) {
        sharedNode node = clique->nodes[i];
        node->colStructure.clear();
        node->colStructure.insert(node->colStructure.end(),
                                  colStructure.begin() + i,
                                  colStructure.end());
    }

    Key lowestReorderedIndex = -1;
    vector<Key> reorderedKeys;
    reorderedKeys.reserve(colStructure.size());

    for(size_t i = 0; i < colStructure.size(); i++) {
        if(lowestReorderedIndex == -1 && colStructure[i] != oldColStructure[i]) {
            lowestReorderedIndex = i;
        }
        if(lowestReorderedIndex != -1) {
            reorderedKeys.push_back(colStructure[i]);
        }
    }

    if(lowestReorderedIndex != -1) {
        // The reordered variables cannot be part of the clique otherwise it would've be reset
        for(sharedNode node : clique->nodes) {
            auto& column = cholesky_.column(node->key) ;

            column.reorderBlocks(reorderedKeys, lowestReorderedIndex);

            // Need to decrement this as the columns in the clique have decreasing number of blocks
            lowestReorderedIndex--;
        }

        for(sharedNode node : clique->nodes) {
            // Reorder factorColStructure only if we need to reorder colStructure
            set<Key, OrderingLess> newFactorColStructure(orderingLess);
            set<Key, OrderingLess> newChangedFactorColStructure(orderingLess);
            for(Key k : node->factorColStructure) {
                newFactorColStructure.insert(k);
            }
            node->factorColStructure = std::move(newFactorColStructure);
            for(Key k : node->changedFactorColStructure) {
                if(!orderingLess(k, node->key)) {
                    newChangedFactorColStructure.insert(k);
                }
            }
            node->changedFactorColStructure = std::move(newChangedFactorColStructure);
        }
    }

    // cout << "reorderedKeys: ";
    // for(const auto& p : reorderedKeys) {
    //     cout << p << " ";
    // }
    // cout << endl;

}

void CholeskyEliminationTree::backwardSolve(VectorValues* delta_ptr, double tol) {
    // cout << "[CholeskyEliminationTree] backwardSolve()" << endl;
    // Do a pre-order traversal from top ot bottom
    // For each node, first process the belowDiagonalBlocks, then do solve on the transpose of the diagonal
    vector<pair<sharedClique, bool>> stack(1, {root_, false});
    while(!stack.empty()) {
        auto& curPair = stack.back();
        sharedClique clique = curPair.first;
        bool& expanded = curPair.second;

        if(!expanded) {
            expanded = true;
            if(clique->orderingVersion() != ordering_version_) {
                reorderClique(clique);
            }

            // cout << "Backsolve clique: ";
            // for(sharedNode node : clique->nodes) {
            //     cout << node->key << " ";
            // }
            // cout << endl;

            bool propagate = false;
            for(Key key : clique->back()->colStructure) {

                if(backSolveKeys_[key]) {
                    propagate = true;
                    break;
                }
            }
            if(propagate) {

                backwardSolveClique(clique, delta_ptr, tol);

                for(sharedClique child : clique->children) {
                    stack.push_back({child, false});
                }
            }
            // else {
            //     cout << "no propagate" << endl;
            // }

        }
        else {
            stack.pop_back();

            for(sharedNode node : clique->nodes) {
                backSolveKeys_[node->key] = false;
            }
        }
    }
    checkInvariant_afterBackSolve();
}

void CholeskyEliminationTree::backwardSolveClique(sharedClique clique, 
    VectorValues* delta_ptr, double tol) {
    // cout << "CholeskyEliminationTree::backwardSolveClique(): ";
    // for(sharedNode node : clique->nodes) {
    //     cout << node->key << " ";
    // }
    // cout << endl;
    // First gather all the columns
    const Key firstKey = clique->front()->key;
    size_t totalWidth = 0;
    size_t totalHeight = 0;
    vector<Key> keys;
    vector<size_t> widths;
    vector<size_t> heights;
    BlockIndexVector blockStartVec;
    vector<LowerTriangularColumnMatrix*> columns;

    ColMajorMatrix m;

    bool gathered = gatherColumns(clique, &m, &totalWidth, &totalHeight, 
                                  &keys, &widths, &heights, &blockStartVec, &columns);
    assert(gathered);

    // Copy over L^-1 Atb row
    Eigen::VectorXd delta = m.block(totalHeight - 1, 0, 1, totalWidth).transpose();
    // delta_.blockRange(firstRow, totalWidth) 
    //     = m.block(totalHeight - 1, 0, 1, totalWidth).transpose();

    if(totalHeight - 1 > totalWidth) {
        // Disregard b row
        const size_t remainingHeight = totalHeight - 1 - totalWidth;
        Eigen::VectorXd gatherX(remainingHeight);
        size_t curRow = 0;
        for(size_t i = 1; i < blockStartVec.size() - 1; i++) {
            // Do not need to gather last row
            size_t otherKey = blockStartVec[i].first;
            size_t otherHeight = blockStartVec[i].second.second;
            gatherX.block(curRow, 0, otherHeight, 1) = delta_.block(otherKey);
            curRow += otherHeight;
        }
        auto B = m.block(totalWidth, 0, remainingHeight, totalWidth); // below diagonal blocks
        // delta_.blockRange(firstRow, totalWidth) -= B.transpose() * gatherX;
        delta -= B.transpose() * gatherX;
        // cout << "B.T = " << endl << B.transpose() << endl;
        // cout << "gathered X = " << gatherX.transpose() << endl;
    }

    // Solve diagonal
    auto D = m.block(0, 0, totalWidth, totalWidth);
    auto LT = D.transpose().triangularView<Eigen::Upper>();
    // LT.solveInPlace(delta_.blockRange(firstRow, totalWidth));
    LT.solveInPlace(delta);

    // scatter results
    size_t curRow = 0;
    for(size_t i = 0; i < keys.size(); i++) {
        // deltaReplaceMask_.push_back(keys[i]);

        for(size_t j = 0; j < widths[i]; j++) {
            double diff = delta_.block(keys[i])(j, 0) - delta(curRow + j, 0);
            // cout << "abs(diff) = " << abs(diff) << endl;
            if(abs(diff) >= tol) {
                // assert(backSolveKeys_[keys[i]] == false);
                backSolveKeys_[keys[i]] = true;
                break;
            }
        }
        delta_.block(keys[i]) = delta.block(curRow, 0, widths[i], 1);
        delta_ptr->at(keys[i]) = delta_.block(keys[i]);
        curRow += widths[i];
    }

    // cout << "delta = " << delta_.blockRange(firstRow, totalWidth).transpose() << endl << endl;

}


/*
void CholeskyEliminationTree::backwardSolveNode(sharedNode node) {
    const Key key = node->key;
    const LowerTriangularColumnMatrix& column = cholesky_.column(key);
    delta_.block(key) = y_.block(key);
    if(column.hasBelowDiagonalBlocks()) {
        const size_t remainingHeight = column.height() - column.width();
        Eigen::VectorXd gatherX(remainingHeight);
        const auto& blockStartVec = column.blockStartVec();
        size_t curRow = 0;
        for(size_t i = 1; i < blockStartVec.size(); i++) {
            size_t otherKey = blockStartVec[i].first;
            size_t otherHeight = blockStartVec[i].second.second;
            gatherX.block(curRow, 0, otherHeight, 1) = delta_.block(otherKey);
            curRow += otherHeight;
        }
        const constBlock B = column.belowDiagonalBlocks();
        delta_.block(key) -= B.transpose() * gatherX;
    }

    // Solve diagonal
    const constBlock D = column.diagonalBlock();
    auto LT = D.transpose().triangularView<Eigen::Upper>();
    LT.solveInPlace(delta_.block(key));

    // cout << "BackSolve. Delta " << key << endl << delta_.block(key) << endl;
}
*/

void CholeskyEliminationTree::updateDelta(VectorValues* delta_ptr) const {
    // cout << "[CholeskyEliminationTree] updateDelta()" << endl;

    for(size_t k = 0; k < delta_ptr->size(); k++) {
        delta_ptr->at(k) = delta_.block(k);
    }
    // delta_ptr->print();
}

bool CholeskyEliminationTree::sorted_no_duplicates(const vector<Key>& v) const {
    if(v.empty()) { return true; }
    Key prev = v.front();
    for(int i = 1; i < v.size(); i++) {
        Key cur = v[i];
        if(!orderingLess(prev, cur)) {
            return false;
        }
        prev = cur;
    }
    return true;
}

void CholeskyEliminationTree::checkInvariant_afterMarkAffected() const {
    for(sharedNode node : nodes_) {
        // assert(sorted_no_duplicates(node->factorColStructure));
        // if(!sorted_no_duplicates(node->changedFactorColStructure)) {
        //     cout << "Node " << node->key << " changedFactorColStructure: ";
        //     for(Key k : node->changedFactorColStructure) {
        //         cout << k << " ";
        //     }
        //     cout << endl;
        // }
        // assert(sorted_no_duplicates(node->changedFactorColStructure));
    }
}

void CholeskyEliminationTree::checkInvariant_afterMarkAncestors() const {
    for(sharedNode node : nodes_) {
        sharedClique clique = node->clique;
        assert(clique != nullptr);
        bool found = false;
        for(sharedNode cliqueNode : clique->nodes) {
            if(node == cliqueNode) {
                found = true;
                break;
            }
        }
        assert(found);
        if(!clique->marked) {
            assert(clique->parent != nullptr);
        }
        else {
            assert(clique->nodes.size() == 1);

            if(clique->parent != nullptr) {
                clique->printClique(cout);
            }
            assert(clique->parent == nullptr);
        }

        for(sharedClique childClique : clique->children) {
            assert(!childClique->marked);
        }
    }
}

void CholeskyEliminationTree::checkInvariant_afterSymbolic() const {
    for(sharedNode node : nodes_) {
        if(node->ordering_version == ordering_version_) {
            assert(sorted_no_duplicates(node->colStructure));
            // assert(sorted_no_duplicates(descendants_.at(node->key)));
            // assert(sorted_no_duplicates(changedDescendants_.at(node->key)));
        }
    }
    assert(root_ != nullptr);
    vector<sharedClique> stack(1, root_);
    vector<bool> reached(nodes_.size(), false);
    while(!stack.empty()) {
        sharedClique clique = stack.back();
        stack.pop_back();

        assert(!clique->nodes.empty());

        vector<Key> nodeKeys;
        for(sharedNode node : clique->nodes) {
            nodeKeys.push_back(node->key);
            assert(reached[node->key] == false);
            reached[node->key] = true;
        }
        assert(sorted_no_duplicates(nodeKeys));

        for(sharedClique childClique : clique->children) {
            stack.push_back(childClique);
        }
    }
    for(bool b : reached) {
        if(!b) {
            for(int i = 0; i < reached.size(); i++) {
                cout << "[" << i << " " << reached[i] << "] ";
            }
            cout << endl;
            exit(0);
        }
        assert(b);
    }

}

void CholeskyEliminationTree::checkInvariant_afterCholesky() const {
    for(auto Status : factorLinearizeStatus_) {
        assert(Status == LINEARIZED);
    }

    for(int i = 0; i < ordering_.size(); i++) {
        assert(is_reordered_[i] == false);
    }

    // validate Tree structure
    vector<bool> key_reached(ordering_.size(), false);
    for(int i = 0; i < ordering_.size(); i++) {
        Key k = orderingToKey_.at(i);
        assert(ordering_[k] == i);
        assert(!key_reached[k]);
        key_reached[k] = true;
    }
    for(int i = 0; i < key_reached.size(); i++) {
        assert(key_reached[i]);
    }

    // Check delta ordering
    for(int i = 0; i < ordering_.size(); i++) {
        auto p = delta_.blockStartVec().at(i);
        Key k = p.first;
        assert(k == orderingToKey_[i]);
        assert(p.second.first == delta_.blockStartMap().at(k).first);
        assert(p.second.second == delta_.blockStartMap().at(k).second);
    }


    // for(const auto& m : descendantStatus) {
    //     for(auto p : m) {
    //         assert(p.second == DONE);
    //     }
    // }

    vector<sharedClique> stack(1, root_);
    vector<bool> node_reached(nodes_.size(), false);
    while(!stack.empty()) {
        sharedClique clique = stack.back();
        stack.pop_back();

        assert(clique->nodes.size() > 0);
        assert(!clique->marked);

        sharedNode lastNode = nullptr;
        for(sharedNode node : clique->nodes) {
            Key key = node->key;
            assert(node_reached[key] == false);
            node_reached[key] = true;

            assert(node->clique == clique);
            if(!node->changedFactorColStructure.empty()) {
                cout << "Node " << node->key << endl;
            }
            assert(node->changedFactorColStructure.empty());
            assert(changedDescendants_.at(node->key).empty());
            assert(!node->is_reordered);
            assert(!node->relinearize);
            assert(*node->colStructure.begin() == node->key);

            // for(Key k : descendants_.at(key)) {
            //     assert(!orderingLess(node->key, k));
            //     bool found = false;
            //     for(Key k : nodes_.at(k)->colStructure) {
            //         if(node->key == k) {
            //             found = true;
            //             break;
            //         }
            //     }
            //     assert(found);
            // }

            auto it = node->factorColStructure.find(node->key);
            while(it != node->factorColStructure.end()) {
                bool found = false;
                for(Key nodeKey : node->colStructure) {
                    if(*it == nodeKey) {
                        found = true;
                        break;
                    }
                }
                assert(found);
                it++;
            }

            if(lastNode != nullptr) {
                // Check clique structure
                assert(lastNode->colStructure.size() == node->colStructure.size() + 1);
                assert(lastNode->colStructure[1] == node->key);
                auto it1 = lastNode->colStructure.begin();
                auto it2 = node->colStructure.begin();
                it1++;
                for(; it1 != lastNode->colStructure.end(); it1++, it2++) {
                    assert(*it1 == *it2);
                }
            }
            lastNode = node;

            // Check cholesky column
            const auto& column = cholesky_.at(node->key);
            assert(column.allocated());
            const auto& blockIndices = column.blockIndices();
            int index = 0;
            for(Key k : node->colStructure) {
                assert(k == blockIndices[index].first);
                index++;
            }
        }

        sharedNode firstNode = clique->nodes.front();
        for(sharedClique childClique : clique->children) {
            assert(childClique->back()->colStructure[1] == firstNode->key);
            assert(childClique->parent == clique);
        }

        for(sharedClique child : clique->children) {
            stack.push_back(child);
        }
    }
    for(int i = 0; i < node_reached.size(); i++) {
        assert(node_reached[i]);
    }
}

void CholeskyEliminationTree::checkInvariant_afterBackSolve() const {
    for(sharedNode node : nodes_) {
        assert(backSolveKeys_[node->key] == false);
    }
}


} // namespace gtsam

// void CholeskyEliminationTree::updateFactorsAndMarkAffectedKeys(
//                                       const NonlinearFactorGraph& nonlinearFactors,
//                                       const FactorIndices& newFactorIndices, 
//                                       const FactorIndices& removeFactorIndices,
//                                       const std::optional<FastList<Key>>& extraKeys,
//                                       KeySet& affectedKeys) {
//     // std::cout << "[CholeskyEliminationTree] updateFactorsAndMarkAffectedKeys()" << std::endl;
//     for(const FactorIndex factorIndex : removeFactorIndices) {
//         sharedFactor factor = nonlinearFactors[factorIndex];
// 
//         for(const Key& key : factor->keys()) {
//             sharedNode node = nodes_[key];
//             node->factorIndexSet.erase(factorIndex);
//             for(const Key& otherKey : factor->keys()) {
//                 assert(node->keyFactorCount[otherKey] > 0);
//                 node->keyFactorCount[otherKey]--;
//             }
//             affectedKeys.insert(key);
//         }
//     }
//     for(const FactorIndex factorIndex : newFactorIndices) {
// 
//         sharedFactor factor = nonlinearFactors[factorIndex];
//         factorLinearizeStatus_[factorIndex] = false;
// 
//         // std::cout << "Factor " << factorIndex << ": ";
//         for(const Key& key : factor->keys()) {
//             assert(key < nodes_.size());
//             // std::cout << key << " ";
//             sharedNode node = nodes_[key];
//             node->factorIndexSet.insert(factorIndex);
// 
//             // Count number of times a node interacts with another node
//             // Only in raw factors
//             for(const Key& otherKey : factor->keys()) {
//                 auto iterPair = node->keyFactorCount.insert({otherKey, 0});
//                 iterPair.first->second++;
//                 hessian_.preallocateOrInitialize(key, otherKey, false);
//             }
//             affectedKeys.insert(key);
// 
//         }
//     }
//     // std::cout << std::endl;
//     if(extraKeys) {
//         for(const Key key : *extraKeys) {
//             affectedKeys.insert(key);
//         }
//     }
//     // exit(1);
//     // TODO:
//     // We need to deal with the smart factor business.
// }
// 
// // TODO: Do not allocate here. Allocate in symbolic elimination
// void CholeskyEliminationTree::markRelinKeys(const KeySet& relinKeys) {
//     // std::cout << "[CholeskyEliminationTree] markRelinKeys()" << std::endl;
//     for(const Key key : relinKeys) {
//         nodes_[key]->relinearized = false;
//         hessian_.resetColumn(key);
//         cholesky_.resetColumn(key);
//     }
// }
// 
// void CholeskyEliminationTree::markAllAffectedKeys(const KeySet& observedKeys, 
//                          const KeySet& relinKeys,
//                          KeySet* markedKeys,
//                          KeyVector* orphanKeys) {
//     // std::cout << "[CholeskyEliminationTree] markAllAffectedKeys()" << std::endl;
//     markedKeys->clear();
//     orphanKeys->clear();
//     for(const Key key : relinKeys) {
//         markNode(nodes_[key], markedKeys);
//         // Mark all keys that have a factor with this key, and mark all their ancestors
//         sharedNode node = nodes_[key];
//         markNode(node);
//         for(const auto& keyCountPair : node->keyFactorCount) {
//             assert(keyCountPair.second > 0);
//             const Key otherKey = keyCountPair.first;
//             markNode(nodes_[otherKey]);
//             // Only reset selected Hessian
//             // Cholesky can be reset in markNode
//             hessian_.preallocateOrInitialize(key, otherKey, true);
//             cholesky_.preallocateOrInitialize(key, otherKey, true);
//             hessian_.preallocateOrInitialize(key, otherKey, true);
//             cholesky_.preallocateOrInitialize(key, otherKey, true);
//         }
//         // Reset corresponding matrices
//     }
// 
//     for(const Key key : observedKeys) {
//         markNode(nodes_[key], markedKeys);
//         // TODO: allocate matrices
//         // All connections with a higher order than this node needs to be 
//         // reset
//         for(const Key otherKey : fillInKeys) {
//             
//         }
//     }
// 
//     // TODO: Mark orphans
//     for(const Key key : *markedKeys) {
//         orphanChildren(nodes_[key], orphanKeys);
//     }
//     // std::cout << "orphan key: ";
//     // for(const Key k : *orphanKeys) {
//     //     std::cout << k << " ";
//     // }
//     // std::cout << std::endl;
// }
// 
// void CholeskyEliminationTree::symbolicElimination(const Ordering& ordering, 
//                                                   const KeyVector& orphanKeys) {
//     // std::cout << "[CholeskyEliminationTree] symbolicElimination()" << std::endl;
//     // ordering.print();
//     // Map from key to ordering for convenience
//     std::unordered_map<Key, size_t> keyToOrdering(ordering.size());
//     size_t i = 0;   // for some reason the ordering class doesn't like random access
//     for(const Key key : ordering) {
//         // TODO: Remove this check. Currently we do this because the set of 
//         // affected keys are different from ISAM
//         if(nodes_[key]->marked) {
//             keyToOrdering.insert({key, i});
//             i++;
//         }
//     }
//     
//     // Need to reparent first to make sure the orphan's fillins get propagated correctly
//     // Orphans reParent
//     for(const Key key : orphanKeys) {
//         sharedNode node = nodes_[key];
//         reParentOrphan(node, keyToOrdering);
//     }
// 
// 
//     size_t myOrder = 0;
//     for(const Key key : ordering) {
//         sharedNode node = nodes_[key];
//         // TODO: Remove this check. Currently we do this because the set of 
//         // affected keys are different from ISAM
//         if(!node->marked) {
// 
//             // assert(node->parent == nullptr);    // Ordering includes child conditionals
//             // std::cout << "Skipping " << key << std::endl;
//             continue;
//         }
//         assert(node->parent == nullptr);
//         assert(node->marked);
//         // std::cout << "sym elim node " << key << std::endl;
// 
//         symbolicEliminateNode(node, myOrder, keyToOrdering);
//         myOrder++;
//     }
// 
//     root_ = nodes_[ordering.back()];   // set the root of the tree to the last node in the ordering
// 
//     validateTree();
// }
// 
// 
// void CholeskyEliminationTree::symbolicEliminateNode(
//         sharedNode node,
//         size_t myOrder,
//         const std::unordered_map<Key, size_t>& keyToOrdering) {
//     node->fillInKeys.clear();
// 
//     // First check all our keyFactorCount for keys that we interact with
//     // Then check our children's fillInKeys
//     // If keys are greater than self in terms of ordering, set own fillInKey
//     for(const auto keyCountPair : node->keyFactorCount) {
//         const Key otherKey = keyCountPair.first;
//         // First check if other key is in ordering
//         // if not, then we can assume that it is our descendant and is not marked
//         auto iter = keyToOrdering.find(otherKey);
//         if(iter == keyToOrdering.end()) {
//             // std::cout << "Not found in ordering. Key = " << otherKey << std::endl;
//             assert(nodes_[otherKey]->marked == false);
//             continue;
//         }
//         // if connected Factor is equal or larger than us, 
//         // there is a nonzero block in the R
//         // matrix. Allocate and initialize
//         size_t otherOrder = iter->second;
//         if(otherOrder < myOrder) { continue; }
// 
//         node->fillInKeys.insert(otherKey);
//     }
// 
//     // now look at all children's fillInKeys
//     for(const sharedNode childNode : node->children) {
//         for(const Key otherKey : childNode->fillInKeys) {
//             // Other key has to be equal to or greater than our key
//             // Other key also has to be in ordering
//             // << " ordering = " << keyToOrdering.at(otherKey) << std::endl;
//             if(otherKey != childNode->key) {
//                 const size_t otherOrder = keyToOrdering.at(otherKey);
//                 assert(otherOrder >= myOrder);
//                 node->fillInKeys.insert(otherKey);
//             }
//         }
//         // TODO: we can use an unordered_map::insert(range) for this
//         // then take out the child key
//     }
// 
//     // Allocate cholesky and contribution matrices
//     // And find parent for node
//     size_t parentOrder = keyToOrdering.size();
//     assert(node->parent == nullptr);
//     for(const Key otherKey : node->fillInKeys) {
//         cholesky_.preallocateOrInitialize(node->key, otherKey, true);
//         // Re-parent here
//         const size_t otherOrder = keyToOrdering.at(otherKey);
//         assert(otherOrder >= myOrder);
//         if(otherOrder > myOrder && otherOrder < parentOrder) {
//             parentOrder = otherOrder;
//             node->parent = nodes_[otherKey];
//         }
//     }
// 
//     if(node->parent) {
//         // connect child to parent
//         node->parent->children.insert(node);
//     }
// 
// }
// 
// void CholeskyEliminationTree::reParentOrphan(
//         sharedNode node, 
//         const std::unordered_map<Key, size_t>& keyToOrdering) {
//     assert(node->parent == nullptr);
//     assert(keyToOrdering.find(node->key) == keyToOrdering.end());
//     size_t parentOrder = keyToOrdering.size();
//     for(const Key otherKey : node->fillInKeys) {
//         if(otherKey == node->key) {
//             continue;
//         }
//         const size_t otherOrder = keyToOrdering.at(otherKey);
//         if(otherOrder < parentOrder) {
//             parentOrder = otherOrder;
//             node->parent = nodes_[otherKey];
//         }
//     }
//     // connect child to parent
//     assert(node->parent != nullptr);
//     node->parent->children.insert(node);
// }
// 
// void CholeskyEliminationTree::resolveAllocate() {
//     hessian_.resolveAllocate();
//     cholesky_.resolveAllocate();
// }
// 
// 
// void CholeskyEliminationTree::choleskyElimination() {
//     // std::cout << "[CholeskyEliminationTree] choleskyElimination()" << std::endl;
//     std::vector<std::pair<sharedNode, bool>> stack(1, {root_, false});
//     while(!stack.empty()) {
//         auto& curPair = stack.back();
//         sharedNode curNode = curPair.first;
//         bool& expanded = curPair.second;
//         if(!expanded) {
//             expanded = true;
//             for(sharedNode child : curNode->children) {
//                 if(child->marked) {
//                     stack.push_back(std::pair<sharedNode, bool>(child, false));
//                 }
//             }
//         }
//         else {
//             choleskyEliminateNode(curNode);
//             curNode->marked = false;
//             // std::cout << "Unmark node " << curNode->key << std::endl;
//             stack.pop_back();
//         }
//     }
// }
// 
// void CholeskyEliminationTree::printTreeStructure(std::ostream& os) {
//     std::vector<sharedNode> stack(1, root_);
//     while(!stack.empty()) {
//         sharedNode curNode = stack.back();
//         stack.pop_back();
//         os << "Node: " << curNode->key << ", Children: ";
//         for(sharedNode child : curNode->children) {
//             os << child->key << " ";
//             stack.push_back(child);
//         }
//         os << std::endl;
//     }
// }
// 
// void CholeskyEliminationTree::markNode(sharedNode node, KeySet* markedKeys) {
//     // std::cout << "Mark node " << node->key;
//     node->marked = true;
//     markedKeys->insert(node->key);
// 
//     sharedNode parent = node->parent;
//     node->parent = nullptr;
// 
//     if(!parent) { return; }
// 
//     // detach node from parent. This will prevent future children nodes from repeating work
//     assert(parent->children.find(node) != parent->children.end());
//     parent->children.erase(node);
// 
//     // Recursively mark parent. This should be tail recursive
//     markNode(parent, markedKeys);
// } 
// 
// 
// void CholeskyEliminationTree::orphanChildren(sharedNode node, KeyVector* orphanKeys) {
//     assert(node->marked); 
//     for(sharedNode child : node->children) {
//         assert(child->parent != nullptr);
//         orphanKeys->push_back(child->key);
//         child->parent = nullptr;
//     }
//     node->children.clear();
// }
// 
// void CholeskyEliminationTree::choleskyEliminateNode(sharedNode node) {
//     // Relinearize       
//     for(const FactorIndex factorIndex : node->factorIndexSet) {
//         
//     }
// }
// 
// void CholeskyEliminationTree::addFillInKey(sharedNode node, const Key otherKey) {
//     const Key key = node->key;
//     sharedNode otherNode = nodes_[otherKey];
// 
//     assert(otherNode->fillInKeys.find(key) == otherNode->fillInKeys.end());
//     assert(node->conditionalKeys.find(key) == node->conditionalKeys.end());
// 
//     node->fillInKeys.insert(otherKey);
//     nodes_[otherKey]->conditionalKeys.insert(node->key);
// }
// 
// void CholeskyEliminationTree::validateTree() {
//     for(sharedNode node : nodes_) {
//         sharedNode parent = node->parent;
//         if(node->parent) {
//             for(Key key : node->fillInKeys) {
//                 if(key == node->key) {
//                     continue;
//                 }
//                 if(parent->fillInKeys.find(key) == parent->fillInKeys.end()) {
//                     std::cout << "Tree is invalid! Node " << node->key << " fill-in " 
//                               << key << " is not in parent " << parent->key << std::endl;
//                     exit(1);
//                 }
//             }
//         }
//     }
// }
// 
// }
