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
public:
    typedef LowerTriangularColumnMatrix::RowHeightPair RowHeightPair;

private:
    CholeskyEliminationTree* eTreePtr = nullptr;

public:

    Key key;
    sharedNode parent = nullptr;
    std::unordered_set<sharedNode> children;

    sharedClique clique = nullptr;

    // 03/21/2023 We're replacing colStructure and colContribution for a simpler design
    // factorColStructure gives us the colStructure from raw factors
    // We use it to determine how many Jacobian columns to load
    // This only include keys that are higher than us in ordering
    // This is because we don't need to load the other A column if it is to the left of us
    // Note: factorColStructure and changedFactorColStructure need to be 
    // reset for reordered variables!
    std::unordered_set<Key> factorColStructure;

    // Columns in A matrix that changed
    std::unordered_set<Key> changedFactorColStructure;
    
    // All colStructure. Need to be sorted so we can determine parent
    // Nodes that are immediate children of reordered nodes also need to run 
    // symbolic factorization to determine new parent
    std::set<Key, OrderingLess> colStructure;

    // Columns that contribute to us
    std::unordered_set<Key> descendants;
    // Columns that contribute to us that got changed
    std::unordered_set<Key> changedDescendants;

    // // Multiset representing the nonzero structure of the column, with the number
    // // representing the number of previous columns
    // // Needs to be a map because we need to know the order
    // std::map<Key, size_t, OrderingLess> colStructure;
    // std::map<Key, size_t, OrderingLess> changedColStructure;

    // // Multiset representing the number of contribution blocks in the column
    // std::unordered_map<Key, size_t> colContribution;
    // // Changed contribution blocks in the column due to relinearization and children blocks
    // // doesn't include new contribution
    // std::unordered_map<Key, size_t> changedColContribution;
    // // New contribution blocks in the column
    // std::unordered_map<Key, size_t> newColContribution;


    // This iteration is reordered
    bool is_reordered = false;

    // This is used to inherit columns
    std::vector<std::pair<Key, RowHeightPair>> lastBlockStartVec;

    // Factors involving this node. Other keys are already stored in the factor
    // Not storing a sharedFactor in case factor got moved
    // This can be a TODO later
    std::vector<FactorIndex> factorIndices;

    // // Columns we're merging as a supernode
    // std::vector<Key> adoptedCols;

    // bool adoptedSkip = false;

    // // Is this column a reconstruct. True by default
    // bool is_reconstruct = true;

    // True if the whole node needs to be relinearized
    bool relinearize = false;

    // // A marked node is one whole entire column will change after eliminimation
    // bool marked = false;

    // Used to check if we need to reorder blocks
    size_t ordering_version = 0;

    Node(CholeskyEliminationTree* eTreePtr_in, const Key key_in)
        : eTreePtr(eTreePtr_in), key(key_in), colStructure(eTreePtr->orderingLess) {
        // colStructure = std::map<Key, size_t, OrderingLess>(eTreePtr);   
        // changedColStructure = std::map<Key, size_t, OrderingLess>(eTreePtr);   
    }

};

} // namespace gtsam
