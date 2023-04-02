/**
* @file    CholeskyEliminationTree.h
* @brief   Elimination tree structure to perform symbolic factorization and Cholesky factorization
* @author  Roger Hsiao (rogerhh)
* @date    Feb. 8, 2023
*/


#include <gtsam/linear/CholeskyEliminationTree.h>
#include <gtsam/linear/CholeskyEliminationTreeNode.h>
#include <gtsam/linear/CholeskyEliminationTreeClique.h>
#include <gtsam/base/LowerTriangularColumnMatrix.h>
#include <gtsam/3rdparty/CCOLAMD/Include/ccolamd.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <algorithm>
#include <fstream>

using namespace std;

namespace gtsam {

using ColMajorMatrix = CholeskyEliminationTree::ColMajorMatrix;

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
        sharedNode newNode = make_shared<Node>(this, key);
        nodes_.push_back(newNode);
        sharedClique newClique = make_shared<Clique>(this);
        newClique->addNode(newNode);
        ordering_.push_back(ordering_.size());
        orderingToKey.push_back(key);
        bool alloc = delta_.preallocateBlock(key, dim, true);
        assert(alloc);
    }
    delta_.resolveAllocate();
} 

void CholeskyEliminationTree::markAffectedKeys(const NonlinearFactorGraph& nonlinearFactors,
                                               const FactorIndices& newFactorIndices,
                                               const KeySet& relinKeys, 
                                               const std::optional<FastList<Key>>& extraKeys,
                                               KeySet* affectedKeys) {
    // cout << "[CholeskyEliminationTree] markAffectedKeys()" << endl;
    affectedKeys->clear();

    // RelinKeys should be processed before we add in factors because we only need to
    // relinearize old factors
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
                    lowerNode->changedFactorColStructure.insert(higherKey);
                }
            }
        }
    }

    for(const FactorIndex factorIndex : newFactorIndices) {
        assert(factorIndex == factorLinearizeStatus_.size());
        factorLinearizeStatus_.push_back(UNLINEARIZED);
        sharedFactor factor = nonlinearFactors[factorIndex];
        factors_.push_back(factor);
        jacobian_.preallocateBlock(factorIndex, -1, factor->dim(), true);

        // Sort factor keys first
        auto factorKeys = factor->keys();
        // Can't sort factor->keys() itself because it messes up linearization
        sort(factorKeys.begin(), factorKeys.end(), orderingLess);
        for(auto it1 = factorKeys.begin(); it1 != factorKeys.end(); it1++) {
            Key lowerKey = *it1;
            sharedNode lowerNode = nodes_[lowerKey];
            lowerNode->factorIndices.push_back(factorIndex);
            affectedKeys->insert(lowerKey);
            jacobian_.preallocateBlock(factorIndex, lowerKey, factor->dim(), true);

            for(auto it2 = it1; it2 != factorKeys.end(); it2++) {
                Key higherKey = *it2;
                // New factors are not part of the changedColStructure
                // As they have to be loaded regardless of edit or reconstruct
                lowerNode->factorColStructure.insert(higherKey);
            }
        }
    }
    jacobian_.resolveAllocate(-1);
}

// Mark all ancestors of directly changed keys, disconnect child from parent and from clique
void CholeskyEliminationTree::markAncestors(const KeySet& affectedKeys, KeySet* markedKeys) {
    // cout << "[CholeskyEliminationTree] markAncestors()" << endl;
    for(const Key key : affectedKeys) {
        markKey(key, markedKeys);
    }
}

void CholeskyEliminationTree::markKey(const Key key, KeySet* markedKeys) {
    if(nodes_[key]->clique->marked) {
        // Node is already processed
        return;
    }

    // cout << "[CholeskyEliminationTree] markKey() " << key << endl;
    sharedNode node = nodes_[key];
    // node->marked = true;
    markedKeys->insert(key);

    // Detach node from clique and orphan that clique
    sharedClique oldClique = node->clique;
    node->clique->detachNode(node);
    // cout << "node " << node->key << " clique->children: ";
    // for(sharedClique clique : node->clique->children) {
    //     cout << clique->back()->key << endl;
    // }
    // cout << endl;

    node->clique->marked = true;
    assert(node->clique != nullptr);
    assert(node->clique->parent == nullptr);


    sharedNode parent = node->parent;
    node->parent = nullptr;  
    if(!parent) { return; }

    assert(parent->children.find(node) != parent->children.end());
    parent->children.erase(node);

    // Recursively mark parent. This should be tail recursive
    markKey(parent->key, markedKeys);
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

    assert(clique->parent == nullptr);
    assert(clique->children.empty());
    assert(clique->nodes.size() == 1);
    assert(clique->marked == false);

    clique->marked = true;

    // Add keys induced by raw factors
    node->colStructure.clear();
    for(const Key higherKey : node->factorColStructure) {
        node->colStructure.insert(higherKey);
    }

    // Merge chilren col structures
    for(sharedNode child : node->children) {
        for(const Key otherKey : child->colStructure) {
            if(otherKey != child->key) {
                node->colStructure.insert(otherKey);
            }
        }
    }

    // Find parent
    if(node->colStructure.size() == 1) {
        assert(ordering_[node->key] == ordering_.size() - 1);
        node->parent = nullptr;
        node->clique->parent = nullptr;
    }
    else {
        auto it = node->colStructure.begin();
        Key parentKey = *(++it);
        node->parent = nodes_[parentKey];
        node->parent->children.insert(node);

        node->clique->setParent(node->parent->clique);
        // cout << "set node " << node->key << " clique parent to node " << node->parent->key << " clique " << node->parent->clique << endl;
    }

    // // DEBUG Checks
    // for(sharedNode child : node->children) {
    //     if(child->parent != node) {
    //         cout << "Child " << child->key << " parent is not node " << node->key << endl;
    //         exit(1);
    //     }
    //     if(child->clique->parent != node->clique) {
    //         cout << "Child " << child->key << " clique parent " << child->clique->parent << " is not node " << node->key << " clique " << node->clique << endl;
    //         exit(1);
    //     }
    // }

    // Check supernode. Only merge if we only have one child
    if(node->children.size() == 1) {
        sharedNode child = *(node->children.begin());
        // We can relax this, can set a number larger than 1
        // But, if we relax this, memory alloc for super column needs to be more careful
        // assert(node->colStructure.size() < child->colStructure.size());
        // Child clique has to also be marked because we're editing or reconstructing
        // cliques together
        if(child->clique->marked 
                && node->colStructure.size() == child->colStructure.size() - 1) {
            // Merge our clique into child clique
            sharedClique oldClique = node->clique;
            mergeCliques(child->clique, node->clique);
            // if(!oldClique.unique()) {
            //     cout << "oldClique should be unique! references: " << oldClique.use_count() << endl;
            //     exit(1);
            // }
            clique = child->clique;
            assert(node->clique == child->clique);
        }
    }

    // node->clique->printClique(cout);

    // Set root
    if(node->colStructure.size() == 1) {
        root_ = node->clique;
    }

    // Set ancestors' descendants
    for(const Key ancestorKey : node->colStructure) {
        sharedNode ancestorNode = nodes_[ancestorKey];
        ancestorNode->descendants.insert(node->key);
        ancestorNode->changedDescendants.insert(node->key);
    }

}

void CholeskyEliminationTree::mergeColStructure(const Key key,
                                                const Key ignoreChildKey,
                                                const map<Key, size_t, OrderingLess>& src, 
                                                map<Key, size_t, OrderingLess>* dest) {
    for(const auto& p : src) {
        const Key otherKey = p.first;
        if(otherKey != ignoreChildKey) {
            assert(key == otherKey || orderingLess(key, otherKey));
            (*dest)[p.first] += p.second;
        }
    }
}

void CholeskyEliminationTree::mergeColContribution(const Key key,
                                                   const Key ignoreChildKey,
                                                   const unordered_map<Key, size_t>& src,
                                                   unordered_map<Key, size_t>* dest) {
    for(const auto& p : src) {
        const Key otherKey = p.first;
        if(otherKey != ignoreChildKey) {
            assert(key == otherKey || orderingLess(key, otherKey));
            (*dest)[p.first] += p.second;
        }
    }
}

void CholeskyEliminationTree::getTotalReordering() {
    // Num of columns
    const size_t nVars = ordering_.size();
    if(nVars == 0) {
        ordering_.clear();
    }
    else if(nVars == 1) {
        ordering_ = vector<size_t>(1, 0);
    }
    // Number of nonzeros in A
    const size_t nEntries = nEntries_;
    // Number of rows
    const size_t nFactors = factorLinearizeStatus_.size();
    const size_t Alen = ccolamd_recommended((int) nEntries, (int) nFactors, (int) nVars);
    vector<int> A = vector<int>(Alen);
    vector<int> p = vector<int>(nVars + 1);
    vector<int> cmember(nVars, 0);
    cmember[nVars - 1] = 1;
    cmember[nVars - 2] = 1;

    p[0] = 0;
    int count = 0;
    KeyVector keys(nVars);
    size_t index = 0;
    // Arrange factors indices into COLAMD format (column major)
    for(sharedNode node : nodes_) {
        const Key key = node->key;
        assert(key == index);
        for(FactorIndex factorIndex : node->factorIndices) {
            A[count++] = (int) factorIndex; // Copy sparse column
        }
        p[++index] = count;   // column j (base 1) goes from A[j-1] to A[j]-1
    }
    assert(count == nEntries);

    //double* knobs = nullptr; /* colamd arg 6: parameters (uses defaults if nullptr) */
    double knobs[CCOLAMD_KNOBS];
    ccolamd_set_defaults(knobs);
    knobs[CCOLAMD_DENSE_ROW] = -1;
    knobs[CCOLAMD_DENSE_COL] = -1;

    int stats[CCOLAMD_STATS]; /* colamd arg 7: colamd output statistics and error codes */

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

    // cout << "New ordering: ";
    for(int i = 0; i < p.size(); i++) {
        if(p[i] == -1) {
            break;
        }
        ordering_[p[i]] = i;
        // cout << p[i] << " ";
    }
    // cout << endl << endl;

    return;
}

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
    cmember->reserve(factorMap.size());
    // Do reorderKeys first
    for(const Key key : reorderKeys) {
        const set<FactorIndex>& s = factorMap.at(key);
        for(const FactorIndex factorIndex : s) {
            A->push_back(allFactors[factorIndex]);
            count++;
        }
        p->push_back(count);
        // Use constraints to fix nodes that are not reordered
        // amongst against nodes that are reordered
        cmember->push_back(1);
    }

    // Constrain the last two variables
    *(cmember->rbegin()) = 2;
    *(cmember->rbegin() + 1) = 2;

    // Then do fixed keys with constraint 0
    for(const Key key : fixedKeys) {
        const set<FactorIndex>& s = factorMap.at(key);
        for(const FactorIndex factorIndex : s) {
            A->push_back(allFactors[factorIndex]);
            count++;
        }
        p->push_back(count);
        // Use constraints to fix nodes that are not reordered
        // amongst against nodes that are reordered
        cmember->push_back(0);
    }


    *nEntries = count;
    *nFactors = allFactors.size();
    *nVars = factorMap.size();
}

void CholeskyEliminationTree::getPartialReordering(const vector<Key>& reorderKeys,
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
    unordered_set<Key> is_reordered;

    constructCSCMatrix(reorderKeys, &nEntries, &nVars, &nFactors, 
                       &A, &p, &cmember, &is_reordered);

    assert(p.size() == nVars + 1);

    const size_t Alen = ccolamd_recommended((int) nEntries, (int) nFactors, (int) nVars);

    A.resize(Alen);

    //double* knobs = nullptr; /* colamd arg 6: parameters (uses defaults if nullptr) */
    double knobs[CCOLAMD_KNOBS];
    ccolamd_set_defaults(knobs);
    knobs[CCOLAMD_DENSE_ROW] = -1;
    knobs[CCOLAMD_DENSE_COL] = -1;

    int stats[CCOLAMD_STATS]; /* colamd arg 7: colamd output statistics and error codes */

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

    partialOrdering->resize(reorderKeys.size());
    size_t index = 0;
    for(int i = 0; i < p.size(); i++) {
        if(p[i] >= reorderKeys.size()) {
            continue;
        }
        else if(p[i] == -1) {
            break;
        }
        const Key key = reorderKeys[p[i]];
        // if(is_reordered.find(key) != is_reordered.end()) {
            // We only care about the keys in reorderKeys
            partialOrdering->at(index++) = key;
        // }
    }
}

void CholeskyEliminationTree::updateOrdering(const KeySet& markedKeys) {
    // cout << "[CholeskyEliminationTree] updateOrdering()" << endl;
    vector<Key> reorderKeys;
    vector<Key> partialOrdering;
    reorderKeys.insert(reorderKeys.begin(), markedKeys.begin(), markedKeys.end());

    size_t nVars = reorderKeys.size();

    getPartialReordering(reorderKeys, &partialOrdering);

    // Adjust ordering_. Shift all fixed keys to the front
    size_t curIndex = 0;
    size_t lowestReorderedIndex = -1;
    vector<Key> deltaReorderKeys;
    for(size_t i = 0; i < orderingToKey.size(); i++) {
        Key key = orderingToKey[i];
        if(markedKeys.find(key) == markedKeys.end()) {
            orderingToKey[curIndex] = key;
            ordering_[key] = curIndex;
            // cout << "Key " << key << " is now position " << curIndex << endl;;

            if(lowestReorderedIndex == -1 && ordering_[key] != curIndex) {
                lowestReorderedIndex = i;
            }
            if(lowestReorderedIndex != -1) {
                deltaReorderKeys.push_back(key);
            }

            curIndex++;
        }
    }
    ordering_version_++;    // update ordering version

    // Update delta_ ordering.
    deltaReorderKeys.insert(deltaReorderKeys.begin(), 
                            partialOrdering.begin(),
                            partialOrdering.end());
    if(lowestReorderedIndex == -1) {
        lowestReorderedIndex = curIndex;
    }
    // cout << "update ordering lowestReorderedIndex = " << lowestReorderedIndex << endl;
    delta_.reorderBlocks(deltaReorderKeys, lowestReorderedIndex);

    // Reset all affected nodes
    for(size_t i = 0; i < partialOrdering.size(); i++) {
        Key key = partialOrdering[i];
        orderingToKey[curIndex + i] = key;
        ordering_[key] = curIndex + i;
        // cout << "Key " << key << " is now position " << curIndex + i << endl;
    }

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

        node->is_reordered = true;

        // Need to reorder factorColStructure and changedFactorColStructure
        // Need to also change descendants 
        unordered_set<Key> newFactorColStructure;
        unordered_set<Key> newChangedFactorColStructure;
        for(const Key otherKey : node->factorColStructure) {
            sharedNode otherNode = nodes_[otherKey];
            if(!orderingLess(otherKey, key)) {
                // otherKey is greater than or equal to us
                newFactorColStructure.insert(otherKey);
            }
            else {
                // otherKey is less than us
                otherNode->factorColStructure.insert(key);
            }
        }
        node->factorColStructure = std::move(newFactorColStructure);

        for(const Key otherKey : node->changedFactorColStructure) {
            sharedNode otherNode = nodes_[otherKey];
            if(!orderingLess(otherKey, key)) {
                newChangedFactorColStructure.insert(otherKey);
            }
            else {
                otherNode->changedFactorColStructure.insert(key);
            }
        }
        node->changedFactorColStructure = std::move(newChangedFactorColStructure);

        // Fix descendants
        for(const Key key : partialOrdering) {
            node->descendants.erase(key);
            node->changedDescendants.erase(key);
        }
    }
}

void CholeskyEliminationTree::reparentOrphanClique(sharedClique clique) {
    // cout << "in reparent orphan" << endl;
    reorderClique(clique);

    // Find a new parent for the last node
    sharedNode node = clique->back();
    // Remove old parent
    node->parent->children.erase(node);
    auto it = node->colStructure.begin();
    Key parentKey = *(++it);
    node->parent = nodes_[parentKey];
    node->parent->children.insert(node);

    // Merge cliques as well
    clique->setParent(node->parent->clique);
}

/*
void CholeskyEliminationTree::updateOrdering(KeySet* affectedKeys) {
    // cout << "[CholeskyEliminationTree] updateOrdering()" << endl;
    getTotalReordering();

    // TODO: Support partial ordering

    // For each variable, reset all relevant data structures and check all factors
    y_.resetBlocks(false);
    delta_.resetBlocks(false);

    map<size_t, Key> orderingToKey;
    for(size_t key = 0; key < ordering_.size(); key++) {
        orderingToKey.insert({ordering_[key], key});
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
        const Key key = orderingToKey[i];

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
        clique->is_reconstruct = true;
        return;
    }

    unordered_set<Key> cliqueFactorColStructure;
    unordered_set<Key> cliqueChangedFactorColStructure;
    unordered_set<Key> cliqueDescendants;
    unordered_set<Key> cliqueChangedDescendants;
    // TODO: optimize for one node cliques
    for(sharedNode node : clique->nodes) {
        for(const Key key : node->factorColStructure) {
            cliqueFactorColStructure.insert(key);
        }
        for(const Key key : node->changedFactorColStructure) {
            cliqueChangedFactorColStructure.insert(key);
        }
        for(const Key key : node->descendants) {
            cliqueDescendants.insert(key);
        }
        for(const Key key : node->changedDescendants) {
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

    // cout << "factor col structure: ";
    // for(const Key key : cliqueFactorColStructure) {
    //     cout << key << " ";
    // }
    // cout << endl;
    // cout << "changed factor col structure: ";
    // for(const Key key : cliqueChangedFactorColStructure) {
    //     cout << key << " ";
    // }
    // cout << endl;
    // cout << "descendants: ";
    // for(const Key key : cliqueDescendants) {
    //     cout << key << " ";
    // }
    // cout << endl;
    // cout << "changed descendants: ";
    // for(const Key key : cliqueChangedDescendants) {
    //     cout << key << " ";
    // }
    // cout << endl;

    size_t totalCols = cliqueFactorColStructure.size() + cliqueDescendants.size();
    size_t changedCols = cliqueChangedFactorColStructure.size() 
                            + cliqueChangedDescendants.size();
    // cout << "totalCols = " << totalCols << " changedCols = " << changedCols << endl;
    if(totalCols == 0 || changedCols >= 0.5 * totalCols) {
        clique->is_reconstruct = true;
    }
    else {
        clique->is_reconstruct = false;
    }
}

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

            // cout << "cholesky here0" << endl;

            // cout << "Restore pass: ";
            // for(sharedNode node : clique->nodes) {
            //     cout << node->key << " ";
            // }
            // cout << endl;

            // Need to reorder fixed nodes before inheriting
            if(!clique->marked) {
                if(clique->orderingVersion() != ordering_version_) {
                    // Reorder each node in the clique
                    reorderClique(clique);
                }
            }

            // cout << "cholesky here1" << endl;
            // Pass down reconstruct columns and edit columns (edit for marked cliques)
            clique->inheritCols();
            // cout << "cholesky here11" << endl;

            // Only expand to children if we are marked or have a reconstruct column
            // It is possible to be marked but all ancestors are edits
            if(clique->marked || !clique->reconstructCols.empty()) {
                if(clique->marked) {
                    assert(clique->orderingVersion() == ordering_version_);
                    // Decide if edit or reconstruct
            // cout << "cholesky here12" << endl;
                    setEditOrReconstruct(clique);
            // cout << "cholesky here13" << endl;

                    // edit from this clique and restore to linear state if edit
                    editAndRestoreFromClique(clique);

            // cout << "cholesky here2" << endl;

                    if(clique->is_reconstruct) {
                        // Need to add to reconstructCols and editCols in reverse order
                        for(auto it = clique->nodes.rbegin(); 
                                 it != clique->nodes.rend(); it++) {
                            sharedNode node = *it;
                            clique->reconstructCols.push_back(node->key);
                            cholesky_.resetColumn(node->key);
                            // cout << "Pushed " << node->key << " to reconstructCols" << endl;
                        }
                    }
                    else {
                        // Need to add to reconstructCols and editCols in reverse order
                        for(auto it = clique->nodes.rbegin(); 
                                 it != clique->nodes.rend(); it++) {
                            sharedNode node = *it;
                            clique->editCols.push_back(node->key);
                            // cout << "Pushed " << node->key << " to editCols" << endl;
                        }
                        // subtract old linearized blocks
                        prepareEditClique(clique);
                    }
            // cout << "cholesky here3" << endl;

                    // Allocate new blocks for edit or reconstruct
                    for(sharedNode node : clique->nodes) {
                        const Key key = node->key;
                        for(const auto& p : node->colStructure) {
                            const Key otherKey = p;
                            assert(!orderingLess(otherKey, key));
                            // If reconstruct, initialize blocks
                            bool alloc = cholesky_.preallocateBlock(
                                                otherKey, key, 
                                                clique->is_reconstruct);
                            // cout << "preallocate " << otherKey << " " << key << " " << alloc << endl;
                        }
                        jacobian_.resolveAllocate(key);
                        cholesky_.resolveAllocate(key);
            // cout << "cholesky here4" << endl;
                    }
                }
                else {
                    reconstructFromClique(clique);
                }
                // only expand to children if we are marked or have a reconstruct
                // It is possible to be marked but all ancestors are edits
                for(sharedClique child : clique->children) {
                    stack.push_back({child, false});
                }
            }
            // cout << "cholesky here5" << endl;

        }
        else {

            // cout << "Eliminate pass: ";
            // for(sharedNode node : clique->nodes) {
            //     cout << node->key << " ";
            // }
            // cout << endl;

            // Eliminate pass
            stack.pop_back();
            // cout << "cholesky here6" << endl;

            // Unmarked clique already reconstructed
            if(clique->marked) {
                // All factors of this clique should be relinearized
                // Do AtA for each node
            // cout << "cholesky here7" << endl;
                prepareEliminateClique(clique, theta);

                // Eliminiate clique
            // cout << "cholesky here8" << endl;
                eliminateClique(clique);

            // cout << "cholesky here9" << endl;
                // Reset node member variables
                for(sharedNode node : clique->nodes) {
                    node->relinearize = false;
                    node->is_reordered = false;
                    node->changedDescendants.clear();
                    node->changedFactorColStructure.clear();
                    // Make a copy of the blockStartVec for inheriting reconstruct keys
                    const LowerTriangularColumnMatrix& column = cholesky_.column(node->key);
                    node->lastBlockStartVec = column.blockStartVec();
                }
            }
            // cout << "cholesky here10" << endl;

            // Reset member variables
            clique->marked = false;
            clique->reconstructCols.clear();
            clique->editCols.clear();
            clique->is_reconstruct = true;

        }

    }
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

    map<size_t, Key> orderingToKey;
    for(int k = 0; k < ordering_.size(); k++) {
        orderingToKey.insert({ordering_[k], k});
    }
    os << "Ordering: " << endl;
    for(const auto p : orderingToKey) {
        os << p.second << " ";
    }
    os << endl << endl;
    os << "Jacobian:" << endl;
    jacobian_.print(os);
    os << "Jacobian end" << endl << endl;
    os << "Cholesky:" << endl;
    cholesky_.print(os);
    os << "Cholesky end" << endl << endl;
    // os << "b:" << endl;
    // b_.print(os);
    // os << "b end" << endl << endl;
    // os << "y:" << endl;
    // y_.print(os);
    // os << "y end" << endl << endl;
    os << "delta:" << endl;
    delta_.print(os);
    os << "delta end" << endl << endl;

    os << "Delta = Vector Values with " << orderingToKey.size() << "elements" << endl;
    for(size_t i = 0; i < orderingToKey.size(); i++) {
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

    // cout << "After eliminate column" << endl;
    // column.print(cout);
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
        vector<pair<Key, RowHeightPair>>* blockStartVec,
        vector<LowerTriangularColumnMatrix*>* columns) {

    // cout << "gatherColumns" << endl;

    // TODO: Optimize this for when clique only has one column
    // At this point, some new columns may not be allocated nor fully allocated
    *totalWidth = 0;
    for(sharedNode node : clique->nodes) {
        auto column_ptr = &cholesky_.column(node->key);
        if(!column_ptr->allocated()) {
            // cout << "skipped " << node->key << endl;
            break;
        }
        keys->push_back(node->key);
        columns->push_back(column_ptr);
        widths->push_back(column_ptr->width());
        heights->push_back(column_ptr->height());
        *totalWidth += widths->back();
    }

    if(columns->empty()) {
        // Generally the final column is in a clique with the second
        // to last column
        return false;
    }

    *totalHeight = columns->front()->height();   // This needs to account for the b row

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

    // Populate blockStartVec for the new column, it can just be diagonal block
    // plus the rest of the last column's blockStartVec
    *blockStartVec = vector<pair<Key, RowHeightPair>>(1, {keys->front(), {0, *totalWidth}});
    const auto& lastBlockStartVec = columns->back()->blockStartVec();
    curRow = *totalWidth;
    for(size_t i = 1; i < lastBlockStartVec.size(); i++) {
        const auto& p = lastBlockStartVec[i];
        const Key key = p.first;
        const size_t height = p.second.second;
        blockStartVec->push_back({key, {curRow, height}});
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
    if(clique->editCols.empty() && clique->is_reconstruct) {
        // If no edit cols 
        return;
    }

    // First gather all the columns
    const Key firstKey = clique->front()->key;
    size_t totalWidth = 0;
    size_t totalHeight = cholesky_.column(firstKey).height();
    vector<Key> keys;
    vector<size_t> widths;
    vector<size_t> heights;
    vector<pair<Key, RowHeightPair>> blockStartVec;
    vector<LowerTriangularColumnMatrix*> columns;

    ColMajorMatrix m;

    // cout << "edit and restore 0" << endl;
    // Here we do need to care if cannot gather
    bool gathered = gatherColumns(clique, &m, &totalWidth, &totalHeight, 
                                  &keys, &widths, &heights, &blockStartVec, &columns);
    // cout << "edit and restore 1" << endl;

    if(!gathered) {
        // cout << "New clique. nothing to Gather" << endl;
        // exit(1);
        return;
    }

    if(!clique->editCols.empty()) {
    // cout << "edit and restore 2" << endl;
        editFromColumn(m, blockStartVec, clique->editCols);
    }
    if(!clique->is_reconstruct) {
    // cout << "edit and restore 3" << endl;
        restoreColumn(m, totalWidth, totalHeight);
    // cout << "edit and restore 4" << endl;
        scatterColumns(m, widths, heights, columns);
    }
    // cout << "edit and restore 5" << endl;
}

void CholeskyEliminationTree::restoreColumn(ColMajorMatrix& m, 
        size_t totalWidth, size_t totalHeight) {

    Block D = m.block(0, 0, totalWidth, totalWidth);
    Block col = m.block(0, 0, totalHeight, totalWidth);
    // Reset D's upper triangular matrix
    Eigen::VectorXd d = D.diagonal();
    D.triangularView<Eigen::Upper>().setZero();
    D.diagonal() = d;
    // We did L^-1 B^T before, now we want to do L*B^T which is B*L^T
    col *= D.transpose();


}

void CholeskyEliminationTree::editFromClique(sharedClique clique) {
    // First gather all the columns
    // then call edit column
    const Key firstKey = clique->front()->key;
    size_t totalWidth = 0;
    size_t totalHeight = cholesky_.column(firstKey).height();
    vector<Key> keys;
    vector<size_t> widths;
    vector<size_t> heights;
    vector<pair<Key, RowHeightPair>> blockStartVec;
    vector<LowerTriangularColumnMatrix*> columns;

    ColMajorMatrix m;

    gatherColumns(clique, &m, &totalWidth, &totalHeight, 
                  &keys, &widths, &heights, &blockStartVec, &columns);
    
    editFromColumn(m, blockStartVec, clique->editCols);

    // If clique is edit, copy the columns back

    if(!clique->is_reconstruct) {
        scatterColumns(m, widths, heights, columns);
    }
    // if(!clique->is_reconstruct) {
    //     size_t curCol = 0, curRow = 0;
    //     for(size_t i = 0; i < keys.size(); i++) {
    //         columns[i]->blockRange(0, heights[i]) = m.block(curRow, curCol, heights[i], widths[i]);
    //         curCol += widths[i];
    //         curRow += widths[i];
    //     }
    // }
}



void CholeskyEliminationTree::editFromColumn(const ColMajorMatrix& m,
    const vector<pair<Key, RowHeightPair>>& blockStartVec,
    const vector<Key>& editCols) {
    assert(!editCols.empty());
    assert(editCols.back().first != LAST_ROW);

    assert(blockStartVec.size() >= 2);

    // Maybe we should have 2 ways of doing this, one where there are not
    // a lot of edit keys, one where there are
    size_t firstRow = blockStartVec[1].second.first;
    size_t bHeight = m.rows() - firstRow;
    size_t bWidth = blockStartVec[0].second.second;             // diagonal height

    auto b = m.block(firstRow, 0, bHeight, bWidth);
    ColMajorMatrix outerProduct(bHeight, bHeight);
    outerProduct.setZero();

    outerProduct.selfadjointView<Eigen::Lower>().rankUpdate(b);
    
    // Edit columns are in reverse order because we add them in from the top of the tree
    auto editIter = editCols.rbegin();
    auto editEnd = editCols.rend();
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
        const auto& destBlockStartMap = destColumn.blockStartMap();
        findDestBlocks(blockStartVec.begin() + i,
                       blockStartVec.end(),
                       destBlockStartMap,
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

        // This should also have handled updating Atb
        assert(blockStartVec.back().first == LAST_ROW);
    }
}

void CholeskyEliminationTree::reconstructFromClique(sharedClique clique) {
    // First gather all the columns
    // then call reconstruct column
    
    const Key firstKey = clique->front()->key;
    size_t totalWidth = 0;
    size_t totalHeight = cholesky_.column(firstKey).height();
    vector<Key> keys;
    vector<size_t> widths;
    vector<size_t> heights;
    vector<pair<Key, RowHeightPair>> blockStartVec;
    vector<LowerTriangularColumnMatrix*> columns;

    ColMajorMatrix m;

    gatherColumns(clique, &m, &totalWidth, &totalHeight, 
                  &keys, &widths, &heights, &blockStartVec, &columns);

    reconstructFromColumn(m, blockStartVec, clique->reconstructCols);

}

void CholeskyEliminationTree::reconstructFromColumn(const ColMajorMatrix& m,
    const std::vector<std::pair<Key, RowHeightPair>>& blockStartVec,
    const std::vector<Key>& reconstructCols) {
    assert(!reconstructCols.empty());
    assert(reconstructCols.back().first != LAST_ROW);

    assert(blockStartVec.size() >= 2);

    // Maybe we should have 2 ways of doing this, one where there are not
    // a lot of reconstruct keys, one where there are
    size_t firstRow = blockStartVec[1].second.first;
    size_t bHeight = m.rows() - firstRow;
    size_t bWidth = blockStartVec[0].second.second;             // diagonal height

    auto b = m.block(firstRow, 0, bHeight, bWidth);
    ColMajorMatrix outerProduct(bHeight, bHeight);
    outerProduct.setZero();

    outerProduct.selfadjointView<Eigen::Lower>().rankUpdate(b);

    // Edit columns are in reverse order because we add them in from the top of the tree
    auto reconstructIter = reconstructCols.rbegin();
    auto reconstructEnd = reconstructCols.rend();
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
        const auto& destBlockStartMap = destColumn.blockStartMap();
        findDestBlocks(blockStartVec.begin() + i,
                       blockStartVec.end(),
                       destBlockStartMap,
                       &destBlockStart);

        // Start from row 0 for reconstruct since we don't reconstruct the diagonal
        size_t curRow = otherRow - bWidth;
        size_t curCol = otherRow - bWidth;
        for(const auto& p : destBlockStart) {
            const size_t destRow = p.first;
            const size_t destHeight = p.second;
            // For reconstruct, we SUBTRACT the block 
            destColumn.blockRange(destRow, destHeight)
                -= outerProduct.block(curRow, curCol, destHeight, otherHeight);
            curRow += destHeight;
        }

        // This should also have handled updating Atb
        assert(blockStartVec.back().first == LAST_ROW);
    }
}

void CholeskyEliminationTree::eliminateClique(sharedClique clique) {
    // cout << "CholeskyEliminationTree::eliminateClique()" << endl;
    // First gather all the columns
    // then call eliminate column
    const Key firstKey = clique->front()->key;
    size_t totalWidth = 0;
    size_t totalHeight = cholesky_.column(firstKey).height();
    vector<Key> keys;
    vector<size_t> widths;
    vector<size_t> heights;
    vector<pair<Key, RowHeightPair>> blockStartVec;
    vector<LowerTriangularColumnMatrix*> columns;

    ColMajorMatrix m;

    gatherColumns(clique, &m, &totalWidth, &totalHeight, 
                  &keys, &widths, &heights, &blockStartVec, &columns);

    eliminateColumn(m, totalWidth, totalHeight, blockStartVec);

    scatterColumns(m, widths, heights, columns);
}

void CholeskyEliminationTree::eliminateColumn(ColMajorMatrix& m,
        const size_t totalWidth,
        const size_t totalHeight,
        const std::vector<std::pair<Key, RowHeightPair>>& blockStartVec) {

    // cout << "Before eliminate column. m = " << endl << m << endl << endl;  

    size_t diagonalWidth = totalWidth;
    Block D = m.block(0, 0, diagonalWidth, diagonalWidth);
    Eigen::LLT<Eigen::Ref<Block>, Eigen::Lower> llt(D);
    auto L = D.triangularView<Eigen::Lower>();

    if(totalHeight > totalWidth) {
        size_t remainingHeight = totalHeight - totalWidth;
        Block B = m.block(totalWidth, 0, remainingHeight, totalWidth);
        L.solveInPlace(B.transpose());

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
            const auto& destBlockStartMap = destColumn.blockStartMap();
            // cout << "eliminate column findDestBlocks otherKey = " << otherKey << endl; 
            findDestBlocks(blockStartVec.begin() + i,
                    blockStartVec.end(),
                    destBlockStartMap,
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
        const std::unordered_map<Key, RowHeightPair>& destBlockStartMap,
        std::vector<RowHeightPair>* destBlockStart) {

    destBlockStart->clear();

    // cout << "destBlockStartMap: ";
    // for(auto p : destBlockStartMap) {
    //     cout << "[" << p.first << "," << p.second.first << "," << p.second.second << "] "; 
    // }
    // cout << endl;

    size_t lastRow = LAST_ROW - 1;    // Use to merge destination blocks
    for(; start != end; start++) {
        const Key destKey = start->first;
        const auto& p = destBlockStartMap.at(destKey);
        size_t destRow = p.first;
        size_t destHeight = p.second;
        // cout << "In findDestBlocks, destKey = " << destKey 
        //      << " destRow = " << destRow << " lastRow = " << lastRow << endl;
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
    // const size_t width = column.width();
    // Block D = column.diagonalBlock();
    // auto L = D.triangularView<Eigen::Lower>();
    // RowMajorMatrix newD(width, width);
    // newD.setZero();
    // newD.triangularView<Eigen::Lower>() = L;

    // if(column.hasBelowDiagonalBlocks()) {
    //     // We did L^-1 * B^T before, now we want to do L * B^T which is B * L^T
    //     Block B = column.belowDiagonalBlocks();
    //     B *= newD.transpose();
    // }
    // D.noalias() = newD * newD.transpose();

    // // At this point, y = L^-1 (Atb - contribution). And we want to get 
    // // (Atb - contribution) back, so want to multiply by L on the left side
    // // On the left side, which is transpose on the right side
    // y_.block(key).transpose() *= newD.transpose();


    // Subtract factors that are relinearized
    for(const FactorIndex factorIndex : node->factorIndices) {
        // Unlinearized factors don't need to be reset
        if(factorLinearizeStatus_[factorIndex] == RELINEARIZED) {
            sharedFactor factor = factors_[factorIndex];
            // Do reverse AtA here, but just for this column
            for(const Key otherKey : factor->keys()) {
                if(!orderingLess(otherKey, key)) {
                    assert(column.blockExists(otherKey));
                    // SUBTRACT AtA block
                    column.block(otherKey) 
                        -= jacobian_.block(factorIndex, otherKey).transpose() * jacobian_.block(factorIndex, key);
                }
            }
            // SUBTRACT Atb
            column.block(LAST_ROW) 
                -= jacobian_.block(factorIndex, -1).transpose() * jacobian_.block(factorIndex, key);
        }
    }
}

void CholeskyEliminationTree::prepareEliminateClique(sharedClique clique, const Values& theta) {
    for(sharedNode node : clique->nodes) {
        prepareEliminateColumn(node, clique->is_reconstruct, theta);
    }
}

void CholeskyEliminationTree::prepareEliminateColumn(sharedNode node, 
        bool is_reconstruct, const Values& theta) {
    // All factors of this node should be relinearized
    const Key key = node->key;
    for(const FactorIndex factorIndex : node->factorIndices) {
        // cout << "factor " << factorIndex << endl;
        sharedFactor factor = factors_[factorIndex];
        if(!is_reconstruct && factorLinearizeStatus_[factorIndex] == LINEARIZED) {
            // Can skip factor if node is edit and factor has already been linearized
            // But cannot skip if the factor is newly linearized, in that
            // case we have to add AtA
            continue;
        }
        if(factorLinearizeStatus_[factorIndex] != LINEARIZED 
                && factorLinearizeStatus_[factorIndex] != NEWLINEARIZED) {
            // cout << "linearize factor " << factorIndex << endl;
            vector<Matrix> A;
            Vector b(factor->dim());
            factor->linearizeToMatrix(theta, &A, &b);
            for(int i = 0; i < factor->size(); i++) {
                const Key factorKey = factor->keys()[i];
                jacobian_.block(factorIndex, factorKey) = A[i];
            }
            jacobian_.block(factorIndex, -1) = b;
            if(factorLinearizeStatus_[factorIndex] == UNLINEARIZED) {
                factorLinearizeStatus_[factorIndex] = NEWLINEARIZED;
            }
            else {
                factorLinearizeStatus_[factorIndex] = LINEARIZED;
            }
        }
        // Do AtA and Atb here, but just for this column
        // If edit, only need to do factors that are un/relienarized
        LowerTriangularColumnMatrix& column = cholesky_.column(key);
        bool lastKeyFlag = true;   // Used to ungrade from NEWLINEARIZED to LINEARIZED
        for(const Key otherKey : factor->keys()) {
            if(!orderingLess(otherKey, key)) {
                assert(column.blockExists(otherKey));
                // TODO: AtA can probably be done at once for each row
                column.block(otherKey) += jacobian_.block(factorIndex, otherKey).transpose() * jacobian_.block(factorIndex, key);
                // cout << "AtA " << otherKey << " " << key << " = " << endl
                //      << jacobian_.block(factorIndex, otherKey).transpose() * jacobian_.block(factorIndex, key) << endl << endl;
                if(otherKey != key) {
                    // If there is another key behind our key
                    lastKeyFlag = false;
                }
            }
        }

        // Do Atb
        column.block(LAST_ROW) 
            += jacobian_.block(factorIndex, -1).transpose() * jacobian_.block(factorIndex, key); 
        // cout << "Atb = "  << jacobian_.block(factorIndex, -1).transpose() * jacobian_.block(factorIndex, key) << endl;
        if(lastKeyFlag) {
            factorLinearizeStatus_[factorIndex] = LINEARIZED;
        }
    }
}

void CholeskyEliminationTree::reorderClique(sharedClique clique) {
    // cout << "in reorder clique" << endl;
    assert(!clique->marked);

    // Find the keys that have been reordered and the lowest reordered key
    const auto& oldColStructure = clique->front()->colStructure; 
    vector<Key> newKeys;
    newKeys.reserve(oldColStructure.size());

    for(const Key key : oldColStructure) {
        newKeys.push_back(key);
    }

    sort(newKeys.begin(), newKeys.end(), orderingLess);

    auto it = oldColStructure.begin();
    Key lowestReorderedIndex = -1;
    vector<Key> reorderedKeys;
    reorderedKeys.reserve(newKeys.size());
    for(size_t i = 0; i < newKeys.size(); i++, it++) {
        if(lowestReorderedIndex == -1 && *it != newKeys[i]) {
            lowestReorderedIndex = i;
        }
        if(lowestReorderedIndex != -1) {
            reorderedKeys.push_back(newKeys[i]);
        }
    }

    // cout << "reorderedKeys: ";
    // for(const auto& p : reorderedKeys) {
    //     cout << p << " ";
    // }
    // cout << endl;

    // The reordered variables cannot be part of the clique otherwise it would've be reset
    for(sharedNode node : clique->nodes) {
        assert(node->ordering_version != ordering_version_);
        node->ordering_version = ordering_version_;
        if(lowestReorderedIndex == -1) {
            continue;
        }
        auto& column = cholesky_.column(node->key) ;

        column.reorderBlocks(reorderedKeys, lowestReorderedIndex);
        set<Key, OrderingLess> newColStructure(orderingLess);
        for(const auto& p : node->colStructure) {
            newColStructure.insert(p);
        }
        node->colStructure = std::move(newColStructure);

        // Need to decrement this as the columns in the clique have decreasing number of blocks
        lowestReorderedIndex--;
    }

    for(sharedNode node : clique->nodes) {
        // Also need to sort lastBlockStartVec
        auto& lastBlockStartVec = node->lastBlockStartVec;
        sort(lastBlockStartVec.begin(), 
                lastBlockStartVec.begin() + lastBlockStartVec.size() - 1, 
                orderingLess);
    }
}

void CholeskyEliminationTree::backwardSolve() {
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

            backwardSolveClique(clique);

            if(!clique->backsolveKeys.empty()) {
                for(sharedClique child : clique->children) {
                    stack.push_back({child, false});
                }
            }
        }
        else {
            stack.pop_back();

            // Reset variables related to backsolve
            clique->backsolveKeys.clear();       
        }
    }
}

void CholeskyEliminationTree::backwardSolveClique(sharedClique clique) {
    // cout << "CholeskyEliminationTree::backwardSolveClique(): ";
    // for(sharedNode node : clique->nodes) {
    //     cout << node->key << " ";
    // }
    // cout << endl;
    // First gather all the columns
    const Key firstKey = clique->front()->key;
    size_t totalWidth = 0;
    size_t totalHeight = cholesky_.column(firstKey).height();
    vector<Key> keys;
    vector<size_t> widths;
    vector<size_t> heights;
    vector<pair<Key, RowHeightPair>> blockStartVec;
    vector<LowerTriangularColumnMatrix*> columns;

    ColMajorMatrix m;

    gatherColumns(clique, &m, &totalWidth, &totalHeight, 
                  &keys, &widths, &heights, &blockStartVec, &columns);

    // Copy over L^-1 Atb row
    const size_t firstRow = delta_.blockStartMap().at(firstKey).first;
    delta_.blockRange(firstRow, totalWidth) 
        = m.block(totalHeight - 1, 0, 1, totalWidth).transpose();

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
        delta_.blockRange(firstRow, totalWidth) -= B.transpose() * gatherX;
        // cout << "B.T = " << endl << B.transpose() << endl;
        // cout << "gathered X = " << gatherX.transpose() << endl;
    }

    // Solve diagonal
    auto D = m.block(0, 0, totalWidth, totalWidth);
    auto LT = D.transpose().triangularView<Eigen::Upper>();
    LT.solveInPlace(delta_.blockRange(firstRow, totalWidth));

    // cout << "delta = " << delta_.blockRange(firstRow, totalWidth).transpose() << endl << endl;

    // Propagate backsolve keys
    clique->inheritBacksolveKeys();
    for(sharedNode node : clique->nodes) {
        clique->backsolveKeys.push_back(node->key);
    }
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

void CholeskyEliminationTree::mergeCliques(sharedClique clique1, sharedClique clique2) {
    // Assert the clique is standalone because we generally only merge with lone nodes
    // during symbolic elimination
    assert(clique2->nodes.size() == 1);   

    if(clique2->nodes.size() != 1) {
        cout << "merged clique is not a lone clique!" << endl;
        exit(1);
    }

    // cout << "Merge to clique " << clique1 << ": [";
    // for(sharedNode node : clique1->nodes) {
    //     cout << node->key << " ";
    // }
    // cout << "] [";
    // for(sharedNode node : clique2->nodes) {
    //     cout << node->key << " ";
    // }
    // cout << "]" << endl;

    clique1->addNode(clique2->front());

    // Added clique can only have 1 node, but may have clique children
    // Need to add clique children to merged clique
    // Need to make sure to not create circular references
    vector<sharedClique> clique2Children;
    clique2Children.insert(clique2Children.begin(), 
                           clique2->children.begin(),
                           clique2->children.end());
    for(sharedClique childClique : clique2Children) {
        if(childClique != clique1) {
            childClique->setParent(clique1);
        }
    }

    // clique1 parent must now be set to clique2 parent
    if(clique2->parent) {
        clique1->setParent(clique2->parent);
    }

    // After reassigning parent, no node shoud point to this clique
    clique2->detachParent();
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
