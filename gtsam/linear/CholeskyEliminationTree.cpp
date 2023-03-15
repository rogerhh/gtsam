/**
* @file    CholeskyEliminationTree.h
* @brief   Elimination tree structure to perform symbolic factorization and Cholesky factorization
* @author  Roger Hsiao (rogerhh)
* @date    Feb. 8, 2023
*/


#include <gtsam/linear/CholeskyEliminationTree.h>
#include <gtsam/linear/CholeskyEliminationTreeNode.h>
#include <gtsam/3rdparty/CCOLAMD/Include/ccolamd.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <algorithm>
#include <fstream>

using namespace std;

namespace gtsam {

CholeskyEliminationTree::CholeskyEliminationTree() 
    : b_(0, 1, false), y_(0, 1, false), delta_(0, 1, false), orderingLess(this) { 
}

void CholeskyEliminationTree::addVariables(const Values& newTheta) {
    // cout << "[CholeskyEliminationTree] addVariables()" << endl;
    for(const auto& keyValPair : newTheta) {
        const Key& key = keyValPair.key;
        const Value& val = keyValPair.value;
        const size_t dim = val.dim();
        assert(key == nodes_.size());
        sharedNode newNode = std::make_shared<Node>(this, key);
        nodes_.push_back(newNode);
        ordering_.push_back(ordering_.size());
        // hessian_.addColumn(key, dim);
        cholesky_.addColumn(key, dim);
        jacobian_.addColumn(key, dim);
        bool alloc = y_.preallocateBlock(key, dim, true);
        assert(alloc);
        alloc = delta_.preallocateBlock(key, dim, true);
        assert(alloc);
    }
    y_.resolveAllocate();
    delta_.resolveAllocate();
} 

void CholeskyEliminationTree::markAffectedKeys(const NonlinearFactorGraph& nonlinearFactors,
                                               const FactorIndices& newFactorIndices,
                                               const KeySet& relinKeys, 
                                               const std::optional<FastList<Key>>& extraKeys,
                                               KeySet* affectedKeys) {
    // cout << "[CholeskyEliminationTree] markAffectedKeys()" << endl;

    // RelinKeys should be processed before we add in factors because we only need to
    // relinearize old factors
    // How do you notify another column that's also being relinearized 
    // without double counting. Sol: Count changed blocks on the factor to factor basis
    // cout << "Relin keys: ";
    for(const Key relinKey : relinKeys) {
        sharedNode relinNode = nodes_[relinKey];
        relinNode->relinearize = true;
        // cout << relinKey << " ";
        for(const FactorIndex factorIndex : relinNode->factorIndices) {
            // Each changed factor will result in n(n+1)/2 changed block
            // Where n is the number of variables
            if(factorLinearizeStatus_[factorIndex] != LINEARIZED) {
                // If linear status of that factor is not linearized
                // Then it's either a new factor or a factor that was already counted
                continue;
            }
            factorLinearizeStatus_[factorIndex] = RELINEARIZED;
            sharedFactor factor = nonlinearFactors[factorIndex];
            // Sort the keys in the factor. This should be fast as there 
            // are not that many keys in a factor
            set<Key, OrderingLess> sortedFactorKeys(factor->keys().begin(),
                                                    factor->keys().end(),
                                                    orderingLess);
            // changedColContribution should only apply to the lowest ordered key
            // It would get propagated up to the ancestor keys
            size_t keyCount = 0;
            const size_t keySize = sortedFactorKeys.size();
            sharedNode lowestNode = nodes_[*sortedFactorKeys.begin()];
            // cout << "[";
            for(auto it = sortedFactorKeys.begin(); 
                     it != sortedFactorKeys.end(); it++, keyCount++) {
                const Key otherKey = *it;
                lowestNode->changedColContribution[otherKey] += keySize - keyCount;
                // We technically don't need to mark the other node here
                // As it will be marked when we mark the ancestors
                // as long as we make should to mark the lowest node
                sharedNode otherNode = nodes_[otherKey];
                // cout << otherKey << " ";
            }
            // cout << "] ";
        }

        for(const auto& p : relinNode->keyFactorCount) {
            assert(p.second > 0);
            const Key otherKey = p.first;
            nodes_[otherKey]->marked = true;
            affectedKeys->insert(otherKey);

        }
    }
    // cout << endl << endl;

    // cout << "Observed Keys: ";
    for(const FactorIndex factorIndex : newFactorIndices) {
        assert(factorIndex == factorLinearizeStatus_.size());
        factorLinearizeStatus_.push_back(UNLINEARIZED);
        sharedFactor factor = nonlinearFactors[factorIndex];
        factors_.push_back(factor);
        nEntries_ += factor->keys().size();

        // cout << "[";
        for(const Key key : factor->keys()) {
            assert(key < nodes_.size());
            sharedNode node = nodes_[key];
            node->marked = true;
            node->factorIndices.push_back(factorIndex);
            affectedKeys->insert(key);
            jacobian_.preallocateBlock(factorIndex, key, factor->dim(), true);
            b_.preallocateBlock(factorIndex, factor->dim(), true);

            // cout << key << " ";

            // Count number of times a node interacts with another node
            // Only in raw factors
            for(const Key otherKey : factor->keys()) {
                auto iterPair = node->keyFactorCount.insert({otherKey, 0});
                iterPair.first->second++;
                if(!orderingLess(otherKey, key)) {
                    // If the other key is in our column, it adds to a contribution block
                    // This includes diagonal block
                    // hessian_.preallocateBlock(key, otherKey, false);
                    node->newColContribution[key]++;
                    if(iterPair.second) {
                        // If we have never seen this block before, then set up the block
                        node->changedColStructure.insert({otherKey, 1});
                    }
                }
            }

        }
        // cout << "] ";
    }
    // cout << endl << endl;

    // cout << "Affected Keys: ";
    // for(const Key k : *affectedKeys) {
    //     cout << k << " ";
    // }
    // cout << endl << endl;
    // Update column ordering here
}

// Mark all ancestors of directly changed keys, disconnect child from parent 
void CholeskyEliminationTree::markAncestors(const KeySet& affectedKeys, KeySet* markedKeys) {
    for(const Key key : affectedKeys) {
        markKey(key, markedKeys);
    }
}

void CholeskyEliminationTree::markKey(const Key key, KeySet* markedKeys) {
    sharedNode node = nodes_[key];
    node->marked = true;
    markedKeys->insert(key);

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
    node->marked = true;

    // Compute current changedColContribution from old colStructure
    // After this point changedColContribution store all contribution blocks
    // resulting from our old column structure
    // Also store a copy of the column column structure in vector form
    vector<Key> oldColStructureVec;
    oldColStructureVec.reserve(node->colStructure.size());
    for(const auto& p : node->colStructure) {
        oldColStructureVec.push_back(p.first);
    }

    // cout << "Key: " << key << " old col structure: ";
    // for(auto& p : node->colStructure) {
    //     cout << "[" << p.first << " " << p.second << "] ";
    // }
    // cout << endl;
    if(!oldColStructureVec.empty()) {
        const size_t numRows = oldColStructureVec.size();
        for(size_t rowIndex = 1; rowIndex < numRows; rowIndex++) {
            const Key otherKey = oldColStructureVec[rowIndex];
            node->changedColContribution[otherKey] += numRows - rowIndex;
            assert(node->changedColContribution[otherKey] > 0);
        }
    }

    // Merge our colStructure with ours and our children's changedColStructure 
    // to get new colStructure
    for(sharedNode child : node->children) {
        mergeColStructure(key, child->key, 
                          child->changedColStructure, 
                          &node->changedColStructure);
    }
    mergeColStructure(key, INT_MAX, node->changedColStructure, &node->colStructure);

    // Check supernode. Only merge if we only have one child
    if(node->children.size() == 1) {
        sharedNode child = *(node->children.begin());
        // We can relax this, can set a number larger than 1
        // But, if we relax this, memory alloc for super column needs to be more careful
        // assert(node->colStructure.size() < child->colStructure.size());
        if(node->colStructure.size() == child->colStructure.size() - 1) {
            // cout << "Key = " << node->key << " Our structure: ";
            // for(auto p : node->colStructure) {
            //     cout << p.first << " ";
            // }
            // cout << endl;
            // cout << "Key = " << child->key << " Child structure: ";
            // for(auto p : child->colStructure) {
            //     cout << p.first << " ";
            // }
            // cout << endl;
            child->adoptedSkip = true;
            // Get child's adopted keys and add child
            swap(node->adoptedCols, child->adoptedCols);
            node->adoptedCols.push_back(child->key);
        }
    }

    // cout << "Key: " << key << " col structure: ";
    // for(auto& p : node->colStructure) {
    //     cout << "[" << p.first << " " << p.second << "] ";
    // }
    // cout << endl;

    // Use new colStructure to compute new contribution blocks
    // i.e. contribution blocks resulting from new column blocks 
    if(!node->colStructure.empty()) {
        auto it = node->colStructure.begin();
        it++;
        const size_t numRows = node->colStructure.size();
        const size_t oldNumRows = oldColStructureVec.size();
        size_t rowIndex = 1, oldIndex = 1;
        for(; it != node->colStructure.end(); it++, rowIndex++, oldIndex++) {
            const Key otherKey = it->first;
            assert(numRows > rowIndex);
            size_t newNumBlocks = numRows - rowIndex;
            size_t oldNumBlocks = 0;
            if(oldIndex < oldNumRows && oldColStructureVec[oldIndex] != otherKey) {
                // If cannot find key in the old column structure, this must be a 
                // completely new block
                oldIndex++;
                oldNumBlocks = 0;
            }
            else {
                // We can find key in the old column structure, compare the difference
                // I.e. how many blocks are below this block in the old column structure vs 
                // the new column structure
                oldNumBlocks = oldNumRows - oldIndex;
            }
            assert(node->newColContribution[otherKey] == 0);
            node->newColContribution[otherKey] = newNumBlocks - oldNumBlocks;
        }
    }
    
    // Merge changed contribution blocks and new contribution blocks from children
    for(sharedNode child : node->children) {
        mergeColContribution(key, child->key, 
                             child->changedColContribution, 
                             &node->changedColContribution);
        mergeColContribution(key, child->key, 
                             child->newColContribution, 
                             &node->newColContribution);
    }

    // // Compare old contribution and new contribution to determine if edit or reconstruct
    // if(node->colContribution.size() > 0) {
    //     auto origIt = node->colContribution.find(key);
    //     assert(origIt != node->colContribution.end());
    //     auto changedIt = node->changedColContribution.find(key);
    //     if(changedIt != node->changedColContribution.end()) {
    //         cout << "Changed blocks: " << changedIt->second 
    //              << " Total blocks: " << origIt->second << endl;
    //         assert(changedIt->second <= origIt->second);
    //     }
    //     else {
    //         cout << "Changed blocks: " << 0
    //              << " Total blocks: " << origIt->second << endl;
    //     }
    // } 

    // Comment: 03/08/2023. We should merge after we determine edit or reconstruct
    // // Merge new contribution blocks to self
    // mergeColContribution(key, INT_MAX, node->newColContribution, &node->colContribution);

    // cout << "Key: " << key << " col contribution: ";
    // for(auto& p : node->colContribution) {
    //     cout << "[" << p.first << " " << p.second << "] ";
    // }
    // cout << endl;
    // cout << "Key: " << key << " new col contribution: ";
    // for(auto& p : node->newColContribution) {
    //     cout << "[" << p.first << " " << p.second << "] ";
    // }
    // cout << endl;
    // cout << "Key: " << key << " changed col contribution: ";
    // for(auto& p : node->changedColContribution) {
    //     cout << "[" << p.first << " " << p.second << "] ";
    // }
    // cout << endl;

    if(node->colStructure.size() == 1) {
        assert(ordering_[node->key] == ordering_.size() - 1);
        node->parent = nullptr;
        root_ = node;
    }
    else {
        // Find parent
        auto it = node->colStructure.begin();
        Key parentKey = (++it)->first;
        node->parent = nodes_[parentKey];
        node->parent->children.insert(node);
        // cout << "Key: " << key << " parent: " << parentKey << endl;
    }

}

// void cleanup() {
//     node->changedColStructure.clear();
// }

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
                
                    // Decide if reconstruct or edit
                    // colContribution can be 0 for a new column
                    // Since we only count old columns
                    // New columns are set to reconstruct but doesn't need to be in reconstructCols
                    // or editCols
                    node->is_reconstruct = true;
                    if(node->colContribution[key] != 0) {
                        if(node->changedColContribution[key] > 0.4 * node->colContribution[key]) {
                        // if(true) {
                            // cout << "node marked. push to reconstruct col" << endl;
                            node->reconstructCols.push_back(key);
                            
                            // Need to also initialize for delta
                            y_.setZero(key);
                        }
                        else {
                            // cout << "node marked. push to edit col" << endl;
                            node->editCols.push_back(key);
                            node->is_reconstruct = false;

                            // need to restore column and subtract the old linearized blocks
                            prepareEditColumn(node);
                        }
                    }
                
                    // Merge new contribution blocks to self. Previously done in symbolic elim
                    mergeColContribution(key, INT_MAX, node->newColContribution, &node->colContribution);
                    
                    // Allocate blocks for edit or reconstruct
                    for(const auto& p : node->colStructure) {
                        const Key otherKey = p.first;
                        assert(!orderingLess(otherKey, key));
                        // If reconstruct, initialize blocks
                        bool alloc = cholesky_.preallocateBlock(otherKey, key, node->is_reconstruct);
                    }

                    jacobian_.resolveAllocate(node->key);
                    cholesky_.resolveAllocate(node->key);
                    b_.resolveAllocate();
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
    os << "b:" << endl;
    b_.print(os);
    os << "b end" << endl << endl;
    os << "y:" << endl;
    y_.print(os);
    os << "y end" << endl << endl;
    os << "delta:" << endl;
    delta_.print(os);
    os << "delta end" << endl << endl;
}

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

void CholeskyEliminationTree::prepareEditColumn(sharedNode node) {
    const Key key = node->key;
    // cout << "[CholeskyEliminationTree] prepareEditColumn. key = " << key << endl;
    SparseColumnBlockMatrix& column = cholesky_.column(key);
    const size_t width = column.width();
    Block D = column.diagonalBlock();
    auto L = D.triangularView<Eigen::Lower>();
    RowMajorMatrix newD(width, width);
    newD.setZero();
    newD.triangularView<Eigen::Lower>() = L;

    if(column.hasBelowDiagonalBlocks()) {
        // We did L^-1 * B^T before, now we want to do L * B^T which is B * L^T
        Block B = column.belowDiagonalBlocks();
        B *= newD.transpose();
    }
    D.noalias() = newD * newD.transpose();

    // At this point, y = L^-1 (Atb - contribution). And we want to get 
    // (Atb - contribution) back, so want to multiply by L on the left side
    // On the left side, which is transpose on the right side
    y_.block(key).transpose() *= newD.transpose();


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
                    column.block(otherKey) -= jacobian_.block(factorIndex, otherKey).transpose() * jacobian_.block(factorIndex, key);
                    y_.block(otherKey) -= jacobian_.block(factorIndex, otherKey).transpose() * b_.block(factorIndex);
                }
            }
        }
    }
}

void CholeskyEliminationTree::backwardSolve() {
    // Do a pre-order traversal from top ot bottom
    // For each node, first process the belowDiagonalBlocks, then do solve on the transpose of the diagonal
    vector<sharedNode> stack(1, root_);
    while(!stack.empty()) {
        sharedNode node = stack.back();
        stack.pop_back();
        backwardSolveNode(node);

        for(sharedNode child : node->children) {
            stack.push_back(child);
        }
    }
}

void CholeskyEliminationTree::backwardSolveNode(sharedNode node) {
    const Key key = node->key;
    const SparseColumnBlockMatrix& column = cholesky_.column(key);
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

void CholeskyEliminationTree::updateDelta(VectorValues* delta_ptr) const {

    for(size_t k = 0; k < delta_ptr->size(); k++) {
        delta_ptr->at(k) = delta_.block(k);
    }
    // delta_ptr->print();
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
