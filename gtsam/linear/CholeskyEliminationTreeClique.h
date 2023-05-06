/**
* @file    CholeskyEliminationTreeClique.h
* @brief   A group of fully connected node in the CholeskyEliminationTree
* @author  Roger Hsiao (rogerhh)
* @date    Feb. 8, 2023
*/

#pragma once

#include <gtsam/linear/CholeskyEliminationTree.h>
#include <gtsam/base/LowerTriangularColumnMatrix.h>
#include <utility>
#include <vector>
#include <unordered_set>
#include <iostream>

namespace gtsam {

// <col_data_ptr, col_data_start, col_row_start, r, c>
using GatherColumnSource = std::tuple<std::shared_ptr<std::vector<double>>, 
                                      size_t, size_t, size_t, size_t>;

class CholeskyEliminationTree::Clique : 
    public std::enable_shared_from_this<CholeskyEliminationTree::Clique> {
private:
    CholeskyEliminationTree* eTreePtr = nullptr;

public:

    using BlockIndexVector = LowerTriangularColumnMatrix::BlockIndexVector;

    sharedClique get_ptr() {
        return shared_from_this();
    }

    // indices is the ordered version of nodes
    std::vector<sharedNode> nodes;

    sharedClique parent = nullptr;

    std::set<sharedClique> children;

    bool is_reconstruct = true;
    bool has_reconstruct = false;
    bool marked = false; 

    size_t workspaceIndex = -1;
    BlockIndexVector blockIndices;

    // The clique splitting issue. Currently, if a node in the middle of a clique is marked
    // We would split it into an unmarked clique and marked cliques
    // However, if the unmarked clique is later affected, and it can be merged back into 
    // the marked cliques, we would need to gather column from different locations

    // col_data is where the matrix lives during and after cholesky
    // But after markAncestors, cliques might be broken up, col_data_ptr points to the 
    // underlying vector the matrix currently lives in. If col_data_start == 0,
    // then this clique owns the matrix
    std::vector<GatherColumnSource> gatherSources;

    // std::shared_ptr<std::vector<double>> col_data_ptr = nullptr;
    // size_t col_data_start = 0;
    // size_t col_row_start = 0;   // if col_data_start is not 0, col_row_start denotes the row
    //                             // the clique actually starts in
    // size_t r = 0, c = 0;

    // How much memory is needed at this clique and all its children
    size_t accumSize = 0;
    // How much memory is needed just for this node
    size_t selfSize = 0;

    Clique(CholeskyEliminationTree* eTreePtr_in);

    // add node to clique
    void addNode(sharedNode node);

    // Mark clique starting from lowest key. Detach all nodes from this clique
    // and add into their own cliques. Children cliques of this clique
    // are kept in this clique. Detach parent and return it
    sharedClique markClique(const Key lowestKey, KeySet* markedKeys);

    // Find new parent clique as the lowest nonzero index in any column of the clique
    // not including the diagonal
    void findParent();

    // Set parent to nullptr
    void detachParent();

    void setParent(sharedClique newParent);

    // Merge otherClique into this clique
    void mergeClique(sharedClique otherClique);

    // Reorder underlying matrix. Return false if nothing needs to be done
    bool reorderColumn(BlockIndexVector& newBlockIndices);

    // After splitting a clique, if there is a reordering, we need to deallocate 
    // all old columns as all the marked variables will be reconstruct
    void deallocateExcessCols();

    sharedNode front();
    sharedNode back();

    Key frontKey();
    Key backKey();

    size_t cliqueSize();

    size_t orderingVersion();

    void printClique(std::ostream& os);


};

} // namespace gtsam

