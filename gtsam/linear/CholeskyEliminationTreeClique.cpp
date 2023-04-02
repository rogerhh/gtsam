/**
* @file    CholeskyEliminationTreeClique.cpp
* @brief   A group of fully connected node in the CholeskyEliminationTree
* @author  Roger Hsiao (rogerhh)
* @date    Mar. 16, 2023
*/

#include <gtsam/linear/CholeskyEliminationTreeClique.h>
#include <gtsam/linear/CholeskyEliminationTreeNode.h>
#include <iostream>
#include <cassert>

namespace gtsam {

using namespace std;

void CholeskyEliminationTree::Clique::addNode(sharedNode node) {
    nodes.push_back(node);
    node->clique = shared_from_this();
}

void CholeskyEliminationTree::Clique::detachNode(sharedNode node) {
    // FIXME: seems like we need to update the ordering of the nodes here
    // Update: We probably don't need to as all reordered nodes must be detached
    // Fixed child nodes will be reassigned to a new parent

    // TODO: Need to figure out why this is not working
    // cout << "In detach node clique children: ";
    // for(sharedClique clique : children) {
    //     cout << clique->back()->key << " ";
    // }
    // cout << endl;
    int i;
    for(i = nodes.size() - 1; i >= 0; i--) {
        if(i > 0) {
            // make new clique for node
            sharedClique newClique = make_shared<Clique>(eTreePtr);
            newClique->addNode(nodes[i]);
            if(nodes[i] == node) {
                // The last detached node's clique should have this clique as a child
                setParent(newClique);
                // resize vector to detach all nodes that have new cliques
                nodes.resize(i);
                break;
            }
        }
        else {
            assert(nodes[i] == node);
            assert(node->clique == get_ptr());
            nodes.resize(1);
            detachParent();
        }
    }
}

void CholeskyEliminationTree::Clique::detachParent() {
    if(parent != nullptr) {
        assert(parent->children.find(get_ptr()) != parent->children.end());
        parent->children.erase(get_ptr());
        parent = nullptr;
    }
}

void CholeskyEliminationTree::Clique::setParent(sharedClique newParent) {
    // Detach from old parent
    detachParent();

    // Set new parent
    parent = newParent;
    assert(parent->children.find(get_ptr()) == parent->children.end());
    parent->children.insert(get_ptr());
}

void CholeskyEliminationTree::Clique::inheritCols() {
    if(parent) {
        // This is the new col structure. We want to use the new col structure so that 
        // things get inherited correctlya. If reordered, there will be no edit cols
        const auto& colStructure = nodes.back()->colStructure;
        auto iter = parent->reconstructCols.begin();
        auto iterEnd = parent->reconstructCols.end();
        // Start from back as editCols and reconstructCols are added to in reverse order
        auto colIter = colStructure.rbegin();
        auto colIterEnd = colStructure.rend();

        while(iter != iterEnd && colIter != colIterEnd) {
            // Iterating through in descending order
            if(eTreePtr->orderingLess(*colIter, *iter)) {
                iter++;
            } 
            else if(eTreePtr->orderingLess(*iter, *colIter)) {
                colIter++;
            }
            else {
                assert(*iter != lastKey);
                reconstructCols.push_back(*iter);
                iter++;
                colIter++;
            }
        }
        // cout << "column: ";
        // for(auto p : colStructure) {
        //     cout << p << " ";
        // }
        // cout << endl;
        // cout << "parent reconstruct: ";
        // for(Key k : parent->reconstructCols) {
        //     cout << k << " ";
        // }
        // cout << endl;
        // cout << "inherited: ";
        // for(Key k : reconstructCols) {
        //     cout << k << " ";
        // }
        // cout << endl;

        if(marked) {
            auto iter = parent->editCols.begin();
            auto iterEnd = parent->editCols.end();
            auto colIter = colStructure.rbegin();
            auto colIterEnd = colStructure.rend();

            // Iterate through in descending order
            while(iter != iterEnd && colIter != colIterEnd) {
                if(eTreePtr->orderingLess(*colIter, *iter)) {
                    iter++;
                } 
                else if(eTreePtr->orderingLess(*iter, *colIter)) {
                    colIter++;
                }
                else {
                    assert(*iter != lastKey);
                    editCols.push_back(*iter);
                    iter++;
                    colIter++;
                }
            }
        }
    }
}

void CholeskyEliminationTree::Clique::inheritBacksolveKeys() {
    if(parent) {
        // const auto& lastColStructure = nodes.back()->colStructure;
        const Key lastKey = nodes.back()->key;
        const auto& blockStartVec = eTreePtr->cholesky_.column(lastKey).blockStartVec();
        auto iter = parent->backsolveKeys.begin();
        auto iterEnd = parent->backsolveKeys.end();
        auto colIter = blockStartVec.begin() + 1;   // Diagonal block does not need to be checked
        auto colIterEnd = blockStartVec.end();

        // cout << "lastKey = " << lastKey << endl;
        // cout << "inherit backsolve blockstartvec: ";
        // for(auto p : blockStartVec) {
        //     cout << "[" << p.first << "," << p.second.first << "," << p.second.second << "] ";
        // }
        // cout << endl;
        // cout << "inherit backsolve parent: ";
        // for(auto p : parent->backsolveKeys) {
        //     cout << p << " ";
        // }
        // cout << endl;

        while(iter != iterEnd && colIter->first != LAST_ROW) {
            // Final block in col is going to be -1, so don't add
            if(eTreePtr->orderingLess(*iter, colIter->first)) {
                iter++;
            } 
            else if(eTreePtr->orderingLess(colIter->first, *iter)) {
                colIter++;
            }
            else {
                assert(*iter != lastKey);
                backsolveKeys.push_back(*iter);
                iter++;
                colIter++;
            }
        }
    }
}

shared_ptr<CholeskyEliminationTree::Node> CholeskyEliminationTree::Clique::front() {
    return nodes.front();
}

shared_ptr<CholeskyEliminationTree::Node> CholeskyEliminationTree::Clique::back() {
    return nodes.back();
}

size_t CholeskyEliminationTree::Clique::orderingVersion() {
    return nodes.front()->ordering_version;
}

void CholeskyEliminationTree::Clique::printClique(ostream& os) {
    os << "Clique: ";
    for(sharedNode node : nodes) {
        os << node->key << " ";
    }
    os << endl;
    os << "   Parent: ";
    for(sharedNode node : parent->nodes) {
        os << node->key << " ";
    }
    os << endl;
    os << "   Children: " << endl;
    for(sharedClique clique : children) {
        os << "   - ";
        for(sharedNode node : clique->nodes) {
            os << node->key << " ";
        }
        os << endl;
    }
}

}   // namespace gtsam
