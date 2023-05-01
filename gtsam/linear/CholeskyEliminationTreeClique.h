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
#include <unordered_set>
#include <iostream>

namespace gtsam {

class CholeskyEliminationTree::Clique : 
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
    bool marked = false; // true;     // New cliques are always marked

    Clique(CholeskyEliminationTree* eTreePtr_in);

    // add node to clique
    void addNode(sharedNode node);

    // detach all nodes after the target node into their own cliques
    // Until the first node. The new cliques do not have parent clique
    // The last node (node) will be the parent to this clique
    // The nodes will keep track of the new cliques
    void detachNode(sharedNode node);

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

    sharedNode front();
    sharedNode back();

    Key frontKey();
    Key backKey();

    size_t orderingVersion();

    void printClique(std::ostream& os);


};

} // namespace gtsam

