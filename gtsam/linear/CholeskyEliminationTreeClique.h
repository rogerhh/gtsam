/**
* @file    CholeskyEliminationTreeClique.h
* @brief   A group of fully connected node in the CholeskyEliminationTree
* @author  Roger Hsiao (rogerhh)
* @date    Feb. 8, 2023
*/

#pragma once

#include <gtsam/linear/CholeskyEliminationTree.h>
#include <utility>
#include <vector>
#include <set>
#include <iostream>

namespace gtsam {

struct CholeskyEliminationTree::Clique : 
    public std::enable_shared_from_this<CholeskyEliminationTree::Clique> {
private:
    CholeskyEliminationTree* eTreePtr = nullptr;

public:

    sharedClique get_ptr() {
        return shared_from_this();
    }

    // indices is the ordered version of nodes
    std::vector<sharedNode> nodes;

    sharedClique parent = nullptr;

    std::set<sharedClique> children;

    bool is_reconstruct = true;
    bool marked = false;

    std::vector<Key> reconstructCols;
    std::vector<Key> editCols;
    std::vector<Key> backsolveKeys;

    Clique(CholeskyEliminationTree* eTreePtr_in) : eTreePtr(eTreePtr_in) {}

    // add node to clique
    void addNode(sharedNode node);

    // detach all nodes after the target node into their own cliques
    // Until the first node. The new cliques do not have parent clique
    // The last node (node) will be the parent to this clique
    // The nodes will keep track of the new cliques
    void detachNode(sharedNode node);

    // Set parent to nullptr
    void detachParent();

    void setParent(sharedClique newParent);

    // Pass down reconstruct columns and edit columns from parent clique
    // edit columns and are passed if clique is marked
    // We only want to inherit columns that we own in the already allocated
    // column. Any new column blocks do not need to edit or reconstruct
    // This is because cliques can only be broken up, so our clique 
    // has to have been together last iteration, which means that everything
    // in this clique has the same column structure
    // At reordering, all affected columns are reset, so no column structure
    void inheritCols();

    // Inherit keys that need to be propagated in backsolve
    void inheritBacksolveKeys();

    sharedNode front();
    sharedNode back();

    size_t orderingVersion();

    void printClique(std::ostream& os);

};



} // namespace gtsam

