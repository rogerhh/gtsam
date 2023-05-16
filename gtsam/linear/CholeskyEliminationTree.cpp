/**
* @file    CholeskyEliminationTree.h
* @brief   Elimination tree structure to perform symbolic factorization and Cholesky factorization
* @author  Roger Hsiao (rogerhh)
* @date    Feb. 8, 2023
*/

// GTSAM TODO: 
// 2. Marginalization
// 3. Ability to add factors to variables in the middle
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

// Convenience function
Eigen::Block<Eigen::Map<ColMajorMatrix>> block(
    Eigen::Map<ColMajorMatrix> m, 
    size_t r, size_t c, size_t h, size_t w) {
    return Eigen::Block<Eigen::Map<ColMajorMatrix>>(m, r, c, h, w);
}

CholeskyEliminationTree::CholeskyEliminationTree() 
    : delta_(0, 1), orderingLess(this) { 
}

void CholeskyEliminationTree::addVariables(const Values& newTheta) {
    // cout << "[CholeskyEliminationTree] addVariables() " << newTheta.size() << endl;
    for(const auto& keyValPair : newTheta) {
        const Key& key = keyValPair.key;
        const Value& val = keyValPair.value;
        const size_t dim = val.dim();
        assert(key == nodes_.size());
        widths_.push_back(dim);
        cholesky_.addColumn(key, dim);
        // jacobian_.addColumn(key, dim);
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

    // Data structure used to remove duplicate keys and sort
    unordered_map<Key, set<Key>> sortedFactorKeys;

    for(const FactorIndex factorIndex : newFactorIndices) {
        assert(factorIndex == factorLinearizeStatus_.size());
        factorLinearizeStatus_.push_back(UNLINEARIZED);
        sharedFactor factor = nonlinearFactors[factorIndex];
        factors_.push_back(factor);
        // jacobian_.preallocateBlock(factorIndex, -1, factor->dim(), true);

        for(Key k1 : factor->keys()) {
            sharedNode node1 = nodes_[k1];
            node1->factorIndices.push_back(factorIndex);
            affectedKeys->insert(k1);
            observedKeys->insert(k1);
            // jacobian_.preallocateBlock(factorIndex, k1, factor->dim(), true);
            for(Key k2 : factor->keys()) {
                node1->factorColStructure.insert(k2);
            }
        }

        size_t height = factor->dim();
        size_t width = 1;
        for(Key k1 : factor->keys()) {
            width += colWidth(k1);
        }
        cachedLinearFactors_.push_back(ColMajorMatrix(height, width));

    }

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
        if(curClique->orderingVersion() != ordering_version_) {
            reorderClique(curClique);
        }

        curClique = curClique->markClique(curKey, markedKeys);
        if(curClique) {
            curKey = curClique->frontKey();
        }
    } while(curClique != nullptr);
}

size_t CholeskyEliminationTree::colWidth(const Key key) const {
    return widths_[key];
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

    if(groupCliqueFlag) {
        // Reorder keys to group cliques together
        groupCliqueFlag = false;

        size_t curIndex = ordering_.size();
        vector<sharedClique> stack(1, root_);
        while(!stack.empty()) {
            sharedClique clique = stack.back();
            stack.pop_back();
            assert(clique->marked);
            assert(is_reordered_[clique->frontKey()]);

            for(auto it = clique->nodes.rbegin(); it != clique->nodes.rend(); it++) {
                sharedNode node = *it;
                curIndex--;
                orderingToKey_[curIndex] = node->key;
                ordering_[node->key] = curIndex;

                assert(is_reordered_[node->key]);
                assert(is_reordered_[orderingToKey_[curIndex]]);
            }

            for(sharedClique childClique : clique->children) {
                if(is_reordered_[childClique->frontKey()]) {
                    stack.push_back(childClique);
                }
            }
        }

        vector<Key> reorderKeys;
        reorderKeys.insert(reorderKeys.end(), 
                           orderingToKey_.begin() + lowestReorderedIndex,
                           orderingToKey_.end());
        delta_.reorderBlocks(reorderKeys, lowestReorderedIndex);
    }

    allocateStack();

    // checkInvariant_afterSymbolic();
}

void CholeskyEliminationTree::symbolicEliminateKey(const Key key) {
    // cout << "[CholeskyEliminationTree] symbolicEliminateKey: " << key << endl;

    sharedNode node = nodes_[key];
    sharedClique clique = node->clique;


    assert(clique->orderingVersion() == ordering_version_);
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

    // FIXME: Seems like there is a case where current clique can have multiple children cliques
    // But only merge to one of them. This is going to cause a problem when splitting cliques
    // Check supernode. Only merge if we only have one child and child is marked
    // This is known as a fundamental clique
    if(clique->children.size() == 1) {
        sharedClique childClique = *(clique->children.begin());
        sharedNode child = childClique->back();
        assert(node->colStructure.size() >= child->colStructure.size() - 1);
        if(childClique->marked 
                && node->colStructure.size() == child->colStructure.size() - 1) {
            childClique->mergeClique(clique);
            // Need to update our pointer to the current cique
            assert(node->clique = childClique);
            clique = childClique;
            
        }
    }

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

    // Merge changed descendants with descendants
    auto& changedDescendants = changedDescendants_.at(node->key);
    // assert(sorted_no_duplicates(changedDescendants));
    if(!changedDescendants.empty()) {
        size_t i = 0;
        auto& descendants = descendants_.at(node->key);
        for(; i < descendants.size(); i++) {
            // Find the index in descendants that is the same as the 
            // first element in changed descendants. If not found,
            // i will be set to the end of descendants
            if(!orderingLess(descendants[i], changedDescendants[0])) {
                break;
            }
            else {
                assert(orderingLess(descendants[i], changedDescendants[0]));
            }
        }
        // At this point we just need to insert changedDescendants into descendants
        descendants.resize(i + changedDescendants.size());
        for(size_t j = 0; j < changedDescendants.size(); j++) {
            descendants[i + j] = changedDescendants[j];
        }
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
        if(p[i] == -1) {
            break;
        }
        Key key = reorderKeys[p[i]];
        if(is_reordered_[key]) {
            partialOrdering->at(index++) = key;
        }
    }

}

void CholeskyEliminationTree::updateOrdering(const KeySet& markedKeys, 
                                             const KeySet& observedKeys) {
    // cout << "[CholeskyEliminationTree] updateOrdering()" << endl;
    groupCliqueFlag = true;

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
    lowestReorderedIndex = -1;
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

    for(int i = 0; i < newOrderingToKey.size(); i++) {
        Key key = newOrderingToKey[i];
        ordering_[key] = i;
    }

    orderingToKey_ = std::move(newOrderingToKey);

    ordering_version_++;    // update ordering version

    if(!delta_.allocated()) {
        cout << "delta not allocated!" << endl;
        exit(1);
    }

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

        assert(node->ordering_version != ordering_version_);
        node->ordering_version = ordering_version_;

        // We had a problem here with cliques not inheriting reconstructCols correctly
        // after reset, since inheritCols requires a sorted version of the 
        // old keys in the new ordering. For edits and for reconstruct (no reordering), 
        // this is not a problem because there is guaranteed to be no reordering, and 
        // blockStartVec provides the old column structure
        // For reconstructs (after reordering) we need the old keys

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
            totalCols += descendants_.at(nodeKey).size();
            changedCols += changedDescendants_.at(nodeKey).size();
        }

        if(totalCols == 0 || changedCols >= 0.9 * totalCols) {
            clique->is_reconstruct = true;
        }
        else {
            clique->is_reconstruct = false;
        }
    }

    if(clique->is_reconstruct) {
        for(sharedNode node : clique->nodes) {
            if(markedStatus_[node->key] != NEW) {
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
    }

}

void CholeskyEliminationTree::checkEditOrReconstructClique(
        sharedClique clique,
        MarkedStatus mode,
        vector<Key>* destCols) {

    const auto& colStructure = clique->back()->colStructure;

    destCols->clear();
    destCols->reserve(colStructure.size());

    for(size_t i = 1; i < colStructure.size(); i++) { 
        // ignore last row
        Key key = colStructure[i];
        if(markedStatus_[key] == mode) {
            destCols->push_back(key);
        }
    }
}

void CholeskyEliminationTree::allocateStack() {
    vector<pair<sharedClique, bool>> stack(1, {root_, false});
    while(!stack.empty()) {
        auto& curPair = stack.back();
        sharedClique clique = curPair.first;
        bool& expanded = curPair.second;

        if(!expanded) {

            expanded = true;

            if(clique->marked) {
                setEditOrReconstruct(clique);
            }

            // FIXME: Checking this seems very inefficient
            for(sharedClique childClique : clique->children) {
                if(childClique->marked) {
                    assert(clique->marked);
                    stack.push_back({childClique, false});
                }
                else {
                    // check if child has reconstructCols
                    auto& colStructure = childClique->back()->colStructure;
                    for(auto it = colStructure.rbegin(); it != colStructure.rend(); it++) {
                        if(markedStatus_[*it] == RECONSTRUCT) {
                            childClique->has_reconstruct = true;
                            stack.push_back({childClique, false});
                            break;
                        }
                    }
                }
            }
        }
        else {
            stack.pop_back();

            // Only count marked children clique and cliques that have reconstruct 
            // for mem allocation
            size_t maxChildSize = 0;
            for(sharedClique childClique : clique->children) {
                if((childClique->marked || childClique->has_reconstruct) 
                        && maxChildSize < childClique->accumSize) {
                    maxChildSize = childClique->accumSize;
                }
            }

            // Count our contribution
            size_t cliqueWidth = 0;
            for(sharedNode node : clique->nodes) {
                cliqueWidth += colWidth(node->key);
            }

            size_t subDiagHeight = 1;   // Include Atb row
            for(size_t i = 1; i < clique->back()->colStructure.size(); i++) {
                Key key = clique->back()->colStructure[i];
                assert(key != -1);
                subDiagHeight += colWidth(key);
            }

            clique->selfSize = (cliqueWidth + subDiagHeight) * (cliqueWidth + subDiagHeight);
            clique->accumSize = clique->selfSize + maxChildSize;
        }
    }

    workspace.allocate(root_->accumSize);
}

void CholeskyEliminationTree::deallocateStack() {
    workspace.deallocate();
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

            // cout << "Restore pass: " << *clique << endl;
            
            if(clique->marked) {
                // check has edit
                vector<Key> editCols;
                checkEditOrReconstructClique(clique, EDIT, &editCols);
                bool gathered = false;
                // if has edit
                if(!editCols.empty()) {
                    // cout << "EditCols: ";
                    // for(Key k : editCols) {
                    //     cout << k << " ";
                    // }
                    // cout << endl;

                    allocateAndGatherClique(clique, true, false);
                    gathered = true;
                    editOrReconstructFromClique(clique, editCols, 1);
                }
                if(!clique->is_reconstruct) {
                    if(!gathered) { 
                        allocateAndGatherClique(clique, true, false);
                        gathered = true;
                    }
                    restoreClique(clique);
                }
                else {
                    if(gathered) {
                        resetCliqueColumns(clique);
                    }
                    else {
                        allocateAndGatherClique(clique, true, true);
                    }
                }
            }
            else {
                // Reorder fixed nodes
                if(clique->orderingVersion() != ordering_version_) {
                    // Reorder each node in the clique
                    reorderClique(clique);
                }

                vector<Key> reconstructCols;
                checkEditOrReconstructClique(clique, RECONSTRUCT, &reconstructCols);

                if(!reconstructCols.empty()) {
                    // cout << "ReconstructCols: ";
                    // for(Key k : reconstructCols) {
                    //     cout << k << " ";
                    // }
                    // cout << endl;
                    allocateAndGatherClique(clique, false, false);
                    editOrReconstructFromClique(clique, reconstructCols, -1);
                }
            }

            // only add child is it's marked or has reconstruct
            // It is possible to be marked but all ancestors are edits
            for(sharedClique child : clique->children) {
                if(child->marked || child->has_reconstruct) {
                    stack.push_back({child, false});
                }
            }

        }
        else {

            // cout << "Eliminate pass: " << *clique << endl;
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

                // Merge with parent
                mergeWorkspaceClique(clique);

                // Scatter workspace back into columns
                scatterClique(clique);

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
            else {
                // Add AtA blocks for reconstruct columns that may be connected to unmarked nodes
                prepareEliminateClique(clique, theta);

                // Merge with parent
                mergeWorkspaceClique(clique);

                // Still need to scatter to clear workspace
                scatterClique(clique, false);
            }

            // cout << "before reset clique" << endl;

            // Reset member variables
            clique->marked = false;
            clique->is_reconstruct = true;
            clique->has_reconstruct = false;
            clique->workspaceIndex = -1;
        }
    }


    // checkInvariant_afterCholesky();
}

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

void CholeskyEliminationTree::allocateAndGatherClique(sharedClique clique, bool allocate, bool reconstruct) {

    // cout << "[CholeskyEliminationTree] allocateAndGatherClique()" << *clique << endl;

    BlockIndexVector& blockIndices = clique->blockIndices;
    
    // TODO: Make blockIndices a permanent thing in a clique
    blockIndices.clear();

    size_t curRow = 0;
    for(Key key : clique->front()->colStructure) {
        size_t width = colWidth(key);

        blockIndices.push_back({key, {curRow, width}});

        curRow += width;
    }
    blockIndices.push_back({LAST_ROW, {curRow, 1}});
    curRow += 1;

    size_t totalHeight = blockIndices.back().second.first + blockIndices.back().second.second;
    const size_t cliqueSize = clique->nodes.size();
    size_t diagWidth = blockIndices[cliqueSize].second.first;
    size_t bHeight = totalHeight - diagWidth;
    assert(totalHeight * totalHeight == clique->selfSize);

    clique->workspaceIndex 
        = workspace.push_matrices({{totalHeight, diagWidth}, {bHeight, bHeight}});

    if(!reconstruct) {
        // Gather column from wherever it is living now
        // col_data_ptr could be from a previously split clique and may come from multiple sources
        Eigen::Map<ColMajorMatrix> newCol = workspace.get_matrices(clique->workspaceIndex)[0];
        size_t curCol = 0;
        for(auto& gatherSource : clique->gatherSources) {
            auto&[col_data_ptr, col_data_start, col_row_start, r, c] = gatherSource;

            assert(col_data_start + r * c <= col_data_ptr->size());
            Eigen::Map<ColMajorMatrix> oldCol(col_data_ptr->data() + col_data_start, r, c);

            assert(r - col_row_start <= newCol.rows());
            assert(curCol + c <= newCol.cols());

            block(newCol, curCol, curCol, r - col_row_start - 1, c) 
                = block(oldCol, col_row_start, 0, r - col_row_start - 1, c);

            // Do last row independently
            block(newCol, totalHeight - 1, curCol, 1, c) 
                = block(oldCol, r - 1, 0, 1, c);

            curCol += c;
        }
        // cout << "gathered \n" << newCol << endl << endl;
    }

    // The column of the clique might be currently living in the space of another clique
    // Release that space now. 
    for(auto& gatherSource : clique->gatherSources) {
        auto&[col_data_ptr, col_data_start, col_row_start, r, c] = gatherSource;
        if(col_data_start != 0) {
            // If we don't own the data
            assert(col_data_ptr->size() >= r * c);
            assert(allocate);
            col_data_ptr->resize(col_data_ptr->size() - r * c);
        }
    }

    if(allocate) {
        clique->gatherSources.clear();

        auto col_data_ptr = make_shared<vector<double>>(totalHeight * diagWidth);

        clique->gatherSources.push_back(make_tuple(col_data_ptr, 0, 0, totalHeight, diagWidth));
    }
    else {
        assert(clique->gatherSources.size() == 1);
        auto&[col_data_ptr, col_data_start, col_row_start, r, c] = clique->gatherSources[0];

        assert(r == totalHeight);
        assert(c == diagWidth);
        assert(col_data_start == 0);
        assert(col_row_start == 0);
    }

}

void CholeskyEliminationTree::scatterClique(sharedClique clique, bool scatter) {
    if(scatter) {
        assert(clique->gatherSources.size() == 1);
        auto&[col_data_ptr, col_data_start, col_row_start, r, c] = clique->gatherSources[0];

        assert(col_data_ptr);

        Eigen::Map<ColMajorMatrix> col(col_data_ptr->data(), r, c);
        col = workspace.get_matrices(clique->workspaceIndex)[0];
    }

    workspace.pop_matrices();
}

void CholeskyEliminationTree::resetCliqueColumns(sharedClique clique) {
    Eigen::Map<ColMajorMatrix> col = workspace.get_matrices(clique->workspaceIndex)[0];
    col.setZero();
}

void CholeskyEliminationTree::restoreClique(
    sharedClique clique) {

    auto col = workspace.get_matrices(clique->workspaceIndex)[0];

    const auto& blockIndices = clique->blockIndices;
    size_t diagWidth = blockIndices[clique->cliqueSize()].second.first;
    auto D = Eigen::Block<Eigen::Map<ColMajorMatrix>>(col, 0, 0, diagWidth, diagWidth);

    D.triangularView<Eigen::StrictlyUpper>().setZero();

    col *= D.transpose();
}


void CholeskyEliminationTree::editOrReconstructFromClique(
    sharedClique clique,
    // Eigen::Map<ColMajorMatrix> m,
    // const Workspace::BlockInfo& blockInfo,
    const vector<Key>& destCols,
    double sign) {

    auto DB = workspace.get_matrices(clique->workspaceIndex)[0];
    auto C = workspace.get_matrices(clique->workspaceIndex)[1];

    const size_t cliqueSize = clique->nodes.size();
    const auto& blockIndices = clique->blockIndices;
    size_t totalHeight = blockIndices.back().second.first + 1;
    size_t diagWidth = blockIndices[cliqueSize].second.first;
    size_t bHeight = totalHeight - diagWidth;
    size_t bWidth = diagWidth;             // diagonal height
    size_t cHeight = bHeight, cWidth = bHeight;

    assert(blockIndices.back().first == -1);
    assert(totalHeight == DB.rows());
    assert(cWidth == totalHeight - diagWidth);

    auto B = Eigen::Block<Eigen::Map<ColMajorMatrix>>(DB, diagWidth, 0, bHeight, bWidth);

    bool processGrouped = true;
    size_t destSize = destCols.size();
    size_t totalSize = blockIndices.size() - 1;
    if(destSize * 1.5 < totalSize) {
        processGrouped = false;
    }

    if(processGrouped) {

        C.selfadjointView<Eigen::Lower>().rankUpdate(B, sign);

        auto destIt = destCols.begin();
        auto destEnd = destCols.end();

        // Set unused columns to be 0 since we did everything at once
        for(size_t i = cliqueSize; i < blockIndices.size(); i++) {
            Key key = blockIndices[i].first;   

            if(destIt == destEnd || *destIt != key) {
                size_t col = blockIndices[i].second.first - diagWidth;   
                size_t width = blockIndices[i].second.second;   

                block(C, 0, col, cHeight, width).setZero();
            }
            else if (destIt != destEnd){
                destIt++;
            }
        }

        assert(destIt == destEnd);
    }
    else {
        auto destIt = destCols.begin();
        auto destEnd = destCols.end();

        for(size_t i = cliqueSize; i < blockIndices.size(); i++) {
            Key key = blockIndices[i].first;   

            if(destIt != destEnd && *destIt == key) {
                destIt++;
                size_t col = blockIndices[i].second.first;   
                size_t cCol = col - diagWidth;   
                size_t width = blockIndices[i].second.second;   
                size_t height = cHeight - cCol;
                block(C, cCol, cCol, height, width) += 
                    sign * block(DB, col, 0, height, diagWidth)
                         * block(DB, col, 0, width, diagWidth).transpose();
            }
        }
    }
}

void CholeskyEliminationTree::eliminateClique(sharedClique clique) {

    auto DB = workspace.get_matrices(clique->workspaceIndex)[0];
    auto C = workspace.get_matrices(clique->workspaceIndex)[1];
    const auto& blockIndices = clique->blockIndices;
    const size_t& cliqueSize = clique->nodes.size();

    size_t totalHeight = blockIndices.back().second.first + 1;
    size_t bWidth = blockIndices[cliqueSize].second.first;
    size_t bHeight = totalHeight - bWidth;

    assert(totalHeight == DB.rows());
    assert(bWidth == DB.cols());

    // // DEBUG
    // ColMajorMatrix m(totalHeight, totalHeight);
    // m.block(0, 0, totalHeight, bWidth) = DB;
    // m.block(bWidth, bWidth, bHeight, bHeight) = C;
    // m.triangularView<Eigen::StrictlyUpper>().setZero();
    // // cout << "Before eliminate. m = \n" << m << endl << endl;
    // DB.triangularView<Eigen::StrictlyUpper>().setZero();
    // cout << "Keys: " << *clique << endl << endl;
    // cout << "Before eliminate. col = \n" << DB << endl << endl;

    Eigen::Block<Eigen::Map<ColMajorMatrix>> D = block(DB, 0, 0, bWidth, bWidth);
    Eigen::LLT<Eigen::Ref<Eigen::Block<Eigen::Map<ColMajorMatrix>>>> llt(D);
    if(llt.info() == Eigen::NumericalIssue) {
        cout << "Diagonal block not positive definite!" << endl;
        clique->printClique(cout);
        exit(1);
    }
    auto L = D.triangularView<Eigen::Lower>();
    auto B = block(DB, bWidth, 0, bHeight, bWidth);

    L.solveInPlace(B.transpose());

    if(bHeight != -1) {
        C.selfadjointView<Eigen::Lower>().rankUpdate(B, -1);
    }

    // // DEBUG
    // m.block(0, 0, totalHeight, bWidth) = DB;
    // m.block(bWidth, bWidth, bHeight, bHeight) = C;
    // m.triangularView<Eigen::StrictlyUpper>().setZero();
    // // cout << "After eliminate. m = \n" << m << endl << endl;
    // cout << "After eliminate. col = \n" << DB << endl << endl;
}

void CholeskyEliminationTree::mergeWorkspaceClique(sharedClique clique) {
    sharedClique parent = clique->parent;
    if(!parent) {
        assert(clique == root_);
        return;
    }

    const size_t childCliqueSize = clique->nodes.size();
    const auto& childBlockIndices = clique->blockIndices;
    auto childC = workspace.get_matrices(clique->workspaceIndex)[1];
    const size_t childDiagWidth = childBlockIndices[childCliqueSize].second.first;

    auto parentDB = workspace.get_matrices(parent->workspaceIndex)[0];
    auto parentC = workspace.get_matrices(parent->workspaceIndex)[1];

    const size_t& parentCliqueSize = parent->nodes.size();
    const auto& parentBlockIndices = parent->blockIndices;
    const size_t parentDiagWidth = parentBlockIndices[parentCliqueSize].second.first;

    assert(childBlockIndices.back().first == -1);
    assert(parentBlockIndices.back().first == -1);

    vector<tuple<size_t, size_t, size_t>> matchedBlockIndices; // (childCol, parentCol, width)

    size_t lastChildCol = -1, lastParentCol = -1;
    size_t j = childCliqueSize;
    for(size_t i = 0; i < parentBlockIndices.size(); i++) {
        const auto& pChild = childBlockIndices.at(j);
        Key childKey = pChild.first;
        size_t childCol = pChild.second.first, childWidth = pChild.second.second;

        assert(childCol >= childDiagWidth);

        const auto& pParent = parentBlockIndices.at(i);
        Key parentKey = pParent.first;
        size_t parentCol = pParent.second.first, parentWidth = pParent.second.second;

        if(parentKey == childKey) {
            assert(childWidth == parentWidth);
            j++;

            if(parentCol != lastParentCol || 
                    (lastParentCol == parentDiagWidth)) {
                // Since we separated DB and C, need to make sure we don't group them together
                matchedBlockIndices.push_back({childCol, parentCol, childWidth});
                lastChildCol = childCol + childWidth;
                lastParentCol = parentCol + parentWidth;
            } 
            else {
                assert(childCol == lastChildCol);
                get<2>(matchedBlockIndices.back()) += childWidth;
            }
        }
    }

    assert(j == childBlockIndices.size());

    // Do clique columns first
    size_t i = 0;
    while(i < matchedBlockIndices.size() && get<1>(matchedBlockIndices[i]) < parentDiagWidth) {
        // copy diagonal
        size_t childCol1, parentCol1, width1;
        std::tie(childCol1, parentCol1, width1) = matchedBlockIndices[i];
        childCol1 -= childDiagWidth;
        // block(parentDB, parentCol1, parentCol1, width1, width1).selfadjointView<Eigen::Lower>()
        //     += block(childC, childCol1, childCol1, width1, width1).selfadjointView<Eigen::Lower>();
        // TODO: Only need to add triangular part
        block(parentDB, parentCol1, parentCol1, width1, width1)
            += block(childC, childCol1, childCol1, width1, width1);

        // subdiagonal
        for(size_t j = i + 1; j < matchedBlockIndices.size(); j++) {
            size_t childCol2, parentCol2, width2;
            std::tie(childCol2, parentCol2, width2) = matchedBlockIndices[j];
            childCol2 -= childDiagWidth;
            block(parentDB, parentCol2, parentCol1, width2, width1)
                += block(childC, childCol2, childCol1, width2, width1);
        }
        i++;
    }

    // Do C matrix
    while(i < matchedBlockIndices.size()) {
        // copy diagonal
        size_t childCol1, parentCol1, width1;
        std::tie(childCol1, parentCol1, width1) = matchedBlockIndices[i];
        childCol1 -= childDiagWidth;
        parentCol1 -= parentDiagWidth;
        block(parentC, parentCol1, parentCol1, width1, width1)
            += block(childC, childCol1, childCol1, width1, width1);
        // block(parentDB, parentCol1, parentCol1, width1, width1)
        //     += block(childDB, childCol1, childCol1, width1, width1);

        // subdiagonal
        for(size_t j = i + 1; j < matchedBlockIndices.size(); j++) {
            size_t childCol2, parentCol2, width2;
            std::tie(childCol2, parentCol2, width2) = matchedBlockIndices[j];
            childCol2 -= childDiagWidth;
            parentCol2 -= parentDiagWidth;
            block(parentC, parentCol2, parentCol1, width2, width1)
                += block(childC, childCol2, childCol1, width2, width1);
        }
        i++;
    }

    // cout << "After merge workspace." << endl;
    // parent->printClique(cout);
    // cout << parentDB << endl << endl;
    // cout << parentC << endl << endl;
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

void CholeskyEliminationTree::prepareEliminateClique(sharedClique clique, const Values& theta) {

    auto DB = workspace.get_matrices(clique->workspaceIndex)[0];
    auto C = workspace.get_matrices(clique->workspaceIndex)[1];
    const auto& blockIndices = clique->blockIndices;
    const size_t& cliqueSize = clique->nodes.size();
    const Key firstCliqueAncestor = blockIndices[cliqueSize].first;
    const size_t diagWidth = blockIndices[cliqueSize].second.first;

    unordered_map<Key, RowHeightPair> keyRowMap;

    for(const auto& p : blockIndices) {
        keyRowMap.insert(p);
    }

    for(sharedNode node : clique->nodes) {

        // Find nodes that are relinearized and the lowest node in the factor
        // must be in our clique
        // We want to subtract all factor-node that is RELINEARIZED-EDIT
        // We DON'T want to add all factor-node that is LINEARIZED-EDIT
        // We want to add all factor-node that is LINEARIZED-RECONSTRUCT
        // We DON'T want to add all factor-node that is LINEARIZED-UNMARKED

        for(const FactorIndex factorIndex : node->factorIndices) {

            sharedFactor factor = factors_[factorIndex];

            // Node must be lowest node in factor
            // We need this to ensure there is an allocated entry in the workspace
            // for this AtA bock
            // FIXME: We shouldn't have to do this check everytime
            bool lowestKeyFlag = true;
            for(const Key factorKey : factor->keys()) {
                if(orderingLess(factorKey, node->key)) {
                    lowestKeyFlag = false;
                    break;
                }
            }
            if(!lowestKeyFlag) {
                continue;
            }

            size_t factorDim = factor->dim();
            ColMajorMatrix& aRow = cachedLinearFactors_[factorIndex];

            vector<tuple<size_t, size_t, size_t>> factorBlockIndices;  // (key, col, width)
            size_t curCol = 0;
            for(Key key : factor->keys()) {
                factorBlockIndices.push_back(make_tuple(key, curCol, colWidth(key)));
                curCol += colWidth(key);
            }
            factorBlockIndices.push_back(make_tuple(-1, curCol, 1));

            // Handle edits here
            if(factorLinearizeStatus_[factorIndex] == RELINEARIZED) {
                for(size_t i = 0; i < factorBlockIndices.size() - 1; i++) {
                    Key lowerKey;
                    size_t col1, width1;
                    tie(lowerKey, col1, width1) = factorBlockIndices[i];

                    if(markedStatus_[lowerKey] != EDIT) {
                        continue;
                    }

                    bool inClique = firstCliqueAncestor == -1 
                        || orderingLess(lowerKey, firstCliqueAncestor);
                    const auto& p1 = keyRowMap.at(lowerKey);
                    const size_t& destCol = p1.first;
                    const size_t& destWidth = p1.second;

                    for(auto& t2 : factorBlockIndices) {
                        Key higherKey;
                        size_t col2, width2;
                        tie(higherKey, col2, width2) = t2;
                        if(higherKey != -1 && orderingLess(higherKey, lowerKey)) {
                            continue;
                        }
                        // cout << "edit AtA " << lowerKey << " " << higherKey << endl;

                        const auto& p2 = keyRowMap.at(higherKey);
                        const size_t& destRow = p2.first;
                        const size_t& destHeight = p2.second;

                        if(inClique) {
                            // SUBTRACT AtA block (lastRow is included)
                            block(DB, destRow, destCol, destHeight, destWidth) 
                                -= aRow.block(0, col2, factorDim, width2).transpose() 
                                    * aRow.block(0, col1, factorDim, width1);
                        }
                        else {
                            // SUBTRACT AtA block (lastRow is included)
                            block(C, destRow - diagWidth, destCol - diagWidth, destHeight, destWidth) 
                                -= aRow.block(0, col2, factorDim, width2).transpose() 
                                    * aRow.block(0, col1, factorDim, width1);
                        }
                    }
                }
            }

            bool factorLinearized = factorLinearizeStatus_[factorIndex] == LINEARIZED;

            // Linearize factor if not already linearized
            if(factorLinearizeStatus_[factorIndex] != LINEARIZED) {
                factorLinearizeStatus_[factorIndex] = LINEARIZED;
                vector<Matrix> A;
                Vector b(factor->dim());
                factor->linearizeToMatrix(theta, &A, &b);

                size_t curCol = 0;
                for(int i = 0; i < factorBlockIndices.size() - 1; i++) {
                    Key key;
                    size_t col;
                    size_t width;
                    tie(key, col, width) = factorBlockIndices[i];

                    assert(A[i].cols() == width);
                    aRow.block(0, col, factorDim, width) = A[i];
                }
                aRow.block(0, get<1>(factorBlockIndices.back()), factorDim, 1) = b;
            }

            for(size_t i = 0; i < factorBlockIndices.size() - 1; i++) {
                Key lowerKey;
                size_t col1, width1;
                tie(lowerKey, col1, width1) = factorBlockIndices[i];

                if(factorLinearized && 
                        (markedStatus_[lowerKey] == EDIT || markedStatus_[lowerKey] == UNMARKED)) {
                    // Skip factor-node that is LINEARIZED-EDIT or LINEARIZE-UNMARKED
                    continue;
                }

                bool inClique = firstCliqueAncestor == -1 
                                    || orderingLess(lowerKey, firstCliqueAncestor);
                const auto& p1 = keyRowMap.at(lowerKey);
                const size_t& destCol = p1.first;
                const size_t& destWidth = p1.second;

                for(auto& t2 : factorBlockIndices) {
                    Key higherKey;
                    size_t col2, width2;
                    tie(higherKey, col2, width2) = t2;
                    if(higherKey != -1 && orderingLess(higherKey, lowerKey)) {
                        continue;
                    }

                    const auto& p2 = keyRowMap.at(higherKey);
                    const size_t& destRow = p2.first;
                    const size_t& destHeight = p2.second;
                
                    if(inClique) {
                        // ADD AtA block (lastRow is included)
                        block(DB, destRow, destCol, destHeight, destWidth) 
                            += aRow.block(0, col2, factorDim, width2).transpose() 
                                * aRow.block(0, col1, factorDim, width1);
                    }
                    else {
                        // ADD AtA block (lastRow is included)
                        block(C, destRow - diagWidth, destCol - diagWidth, destHeight, destWidth) 
                            += aRow.block(0, col2, factorDim, width2).transpose() 
                                * aRow.block(0, col1, factorDim, width1);
                    }
                }
            }
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

    BlockIndexVector newBlockIndices;
    newBlockIndices.reserve(clique->blockIndices.size());
    size_t curRow = 0;
    for(Key key : colStructure) {
        size_t height = colWidth(key);
        newBlockIndices.push_back({key, {curRow, height}});
        curRow += height;
    }
    newBlockIndices.push_back({-1, {curRow, 1}});
    curRow += 1;

    bool reordered = clique->reorderColumn(newBlockIndices);

    if(reordered) { 

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
        }
        else {
            stack.pop_back();

            for(sharedNode node : clique->nodes) {
                backSolveKeys_[node->key] = false;
            }
        }
    }

    deallocateStack();

    // checkInvariant_afterBackSolve();
}

void CholeskyEliminationTree::backwardSolveClique(sharedClique clique, 
    VectorValues* delta_ptr, double tol) {
    // cout << "CholeskyEliminationTree::backwardSolveClique(): " << *clique << endl;

    // Columns are already gathered
    assert(clique->gatherSources.size() == 1);
    auto&[col_data_ptr, col_data_start, col_row_start, r, c] = clique->gatherSources[0];

    assert(col_data_ptr);
    assert(col_data_ptr->size() == r * c);
    
    const auto& blockIndices = clique->blockIndices;
    size_t cliqueSize = clique->nodes.size();
    size_t diagWidth = c;
    size_t totalHeight = r;
    assert(diagWidth == blockIndices[cliqueSize].second.first);

    Eigen::Map<ColMajorMatrix> col(col_data_ptr->data(), r, c);

    // cout << "backsolve gathered col = \n" << col << "\n\n"; 

    // Copy over L^-1 Atb row
    size_t firstRow = delta_.blockStartMap().at(clique->frontKey()).first;
    Eigen::VectorXd oldDelta = delta_.blockRange(firstRow, diagWidth);
    Block delta = delta_.blockRange(firstRow, diagWidth);
    delta = block(col, totalHeight - 1, 0, 1, diagWidth).transpose();

    if(totalHeight - 1 > diagWidth) {
        // Disregard b row
        const size_t remainingHeight = totalHeight - 1 - diagWidth;
        Eigen::VectorXd gatherX(remainingHeight);
        size_t curRow = 0;
        for(size_t i = cliqueSize; i < blockIndices.size() - 1; i++) {
            // Do not need to gather last row
            size_t otherKey = blockIndices[i].first;
            size_t otherHeight = blockIndices[i].second.second;
            gatherX.block(curRow, 0, otherHeight, 1) = delta_.block(otherKey);
            curRow += otherHeight;
        }

        auto B = block(col, diagWidth, 0, remainingHeight, diagWidth); // sub-diagonal blocks
        delta -= B.transpose() * gatherX;
        // cout << "B.T = " << endl << B.transpose() << endl;
        // cout << "gathered X = " << gatherX.transpose() << endl;
    }

    // Solve diagonal
    auto D = block(col, 0, 0, diagWidth, diagWidth);
    auto LT = D.transpose().triangularView<Eigen::Upper>();
    LT.solveInPlace(delta);

    size_t curRow = 0;
    for(sharedNode node : clique->nodes) {
        Key key = node->key;
        size_t width = colWidth(key);
        delta_ptr->at(key) = delta.block(curRow, 0, width, 1);
        for(size_t j = 0; j < width; j++) {
            double diff = delta(curRow + j, 0) - oldDelta(curRow + j, 0);
            // cout << "abs(diff) = " << abs(diff) << endl;
            if(abs(diff) >= tol) {
                backSolveKeys_[key] = true;
                break;
            }
        }
        curRow += width;
    }

}

void CholeskyEliminationTree::updateDelta(VectorValues* delta_ptr) const {
    // cout << "[CholeskyEliminationTree] updateDelta()" << endl;

    for(size_t k = 0; k < delta_ptr->size(); k++) {
        delta_ptr->at(k) = delta_.block(k);
    }
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
    for(size_t i = 0; i < ordering_.size(); i++) {
        Key key = orderingToKey_[i];
        if(i > 0) {
            assert(orderingLess(orderingToKey_[i - 1], orderingToKey_[i]));
        }
        assert(i == ordering_[key]);
    }

    for(sharedNode node : nodes_) {
        if(node->ordering_version == ordering_version_) {
            assert(sorted_no_duplicates(node->colStructure));
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
