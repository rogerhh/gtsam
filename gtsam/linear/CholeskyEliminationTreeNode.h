/**
* @file    CholeskyEliminationTreeNode.h
* @brief   A node in the CholeskyEliminationTree
* @author  Roger Hsiao (rogerhh)
* @date    Feb. 8, 2023
*/

#pragma once

#include <gtsam/linear/CholeskyEliminationTree.h>
#include <utility>

namespace gtsam {

class CholeskyEliminationTree::Node {
private:
    CholeskyEliminationTree* eTreePtr = nullptr;

public:

    Key key;
    sharedNode parent = nullptr;
    std::unordered_set<sharedNode> children;

    // Multiset representing the nonzero structure of the column, with the number
    // representing the number of previous columns
    // Needs to be a map because we need to know the order
    std::map<Key, size_t, OrderingLess> colStructure;
    std::map<Key, size_t, OrderingLess> changedColStructure;

    // Multiset representing the number of contribution blocks in the column
    std::unordered_map<Key, size_t> colContribution;
    // Changed contribution blocks in the column due to relinearization and children blocks
    // doesn't include new contribution
    std::unordered_map<Key, size_t> changedColContribution;
    // New contribution blocks in the column
    std::unordered_map<Key, size_t> newColContribution;

    // Number of factors we share we any particular key
    // Need to keep track of number in case of a factor removal
    std::unordered_map<Key, size_t> keyFactorCount;

    // Factors involving this node. Other keys are already stored in the factor
    // Not storing a sharedFactor in case factor got moved
    // This can be a TODO later
    std::vector<FactorIndex> factorIndices;

    // Columns we're merging as a supernode
    std::vector<Key> adoptedCols;

    bool adoptedSkip = false;

    // Is this column a reconstruct. True by default
    bool is_reconstruct = true;

    // True if the whole node needs to be relinearized
    bool relinearize = false;

    // A marked node is one whole entire column will change after eliminimation
    bool marked = false;

    // A force_reconstruct is determined at symbolic elim time (used for reordering)
    // After force_reconstruct, all ancestors of the node are assumed to be reconstruct
    // Do not check edit, do not need to pass down edit map or reconstruct map
    bool forceReconstruct = false;

    // A batch_reconstruct is determined at restore pass. After batch_reconstruct
    // all children are assumed to be reconstruct. Do not need to check self edit
    // Only pass down edit map to check which ancestors need to be edited
    // Hopefully this edit map will be small
    bool batchReconstruct = false;

    // allReconstruct guarantees that all ancestors of this node is a reconstruct
    // Thus eliminating the need to check for edits
    bool allReconstruct = false;

    // allEdit guarantees that all ancestors of this node is an edit, you can 
    bool allEdit = false;

    // If there are new column blocks that are not placed at the end
    bool reorderNeeded = false;

    std::vector<Key> reconstructCols;
    std::vector<Key> editCols;

    Node(CholeskyEliminationTree* eTreePtr_in, const Key key_in) 
        : eTreePtr(eTreePtr_in), key(key_in) {
        colStructure = std::map<Key, size_t, OrderingLess>(eTreePtr);   
        changedColStructure = std::map<Key, size_t, OrderingLess>(eTreePtr);   
    }

    // Done with elim, reset all appropriate variables
    void done();

};

} // namespace gtsam
