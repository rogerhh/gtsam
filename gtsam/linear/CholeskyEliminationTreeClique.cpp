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

using namespace std;

namespace gtsam {

using sharedNode = CholeskyEliminationTree::sharedNode;
using sharedClique = CholeskyEliminationTree::sharedClique;

CholeskyEliminationTree::Clique::Clique(CholeskyEliminationTree* eTreePtr_in)
    : eTreePtr(eTreePtr_in) {}

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

sharedClique CholeskyEliminationTree::Clique::markClique(const Key lowestKey, KeySet* markedKeys) {
    // cout << "markClique ";
    // for(sharedNode node : nodes) {
    //     cout << node->key << " ";
    // }
    // cout << endl;
    sharedClique oldParent = parent;
    detachParent();
    int i;
    for(i = nodes.size() - 1; i >= 0; i--) {
        markedKeys->insert(nodes[i]->key);
        if(i > 0) {
            // make new clique for node
            sharedClique newClique = make_shared<Clique>(eTreePtr);
            newClique->addNode(nodes[i]);
            newClique->marked = true;
            if(nodes[i]->key == lowestKey) {
                // The last detached node's clique should have this clique as a child
                setParent(newClique);
                // resize vector to detach all nodes that have new cliques
                nodes.resize(i);
                break;
            }
        }
        else {
            marked = true;
            assert(nodes[i]->key == lowestKey);
            assert(nodes[i]->clique == get_ptr());
            nodes.resize(1);
        }
    }
    return oldParent;
}

void CholeskyEliminationTree::Clique::findParent() {
    const auto& colStructure = nodes.back()->colStructure;
    if(colStructure.size() > 1) {
        Key parentKey = colStructure[1];
        assert(eTreePtr->nodes_[parentKey]->clique != nullptr);
        setParent(eTreePtr->nodes_[parentKey]->clique);
        assert(parent->children.size() > 0);
    }
    else {
        parent = nullptr;
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
    if(parent != nullptr) {
        assert(parent->children.find(get_ptr()) == parent->children.end());
        parent->children.insert(get_ptr());
    }
}

shared_ptr<CholeskyEliminationTree::Node> CholeskyEliminationTree::Clique::front() {
    return nodes.front();
}

shared_ptr<CholeskyEliminationTree::Node> CholeskyEliminationTree::Clique::back() {
    return nodes.back();
}

Key CholeskyEliminationTree::Clique::frontKey() {
    return nodes.front()->key;
}

Key CholeskyEliminationTree::Clique::backKey() {
    return nodes.back()->key;
}

size_t CholeskyEliminationTree::Clique::orderingVersion() {
    size_t ordering_version = nodes.front()->ordering_version;
    for(sharedNode node : nodes) {
        if(node->ordering_version != ordering_version) {
            cout << "ordering version in clique mismatched! First node version: " << ordering_version << endl;
            for(sharedNode node : nodes) {
                cout << "node: " << node->key << " version: " << node->ordering_version << endl;
            }
        }
        assert(node->ordering_version == ordering_version);
    }
    return nodes.front()->ordering_version;
}

void CholeskyEliminationTree::Clique::mergeClique(sharedClique otherClique) {
    // Assert the clique is standalone because we generally only merge with lone nodes
    // during symbolic elimination
    assert(parent = otherClique);
    assert(otherClique->nodes.size() == 1);   
    assert(otherClique->children.find(get_ptr()) != otherClique->children.end());   

    addNode(otherClique->front());

    // Added clique can only have 1 node, but may have clique children
    // Need to add clique children to merged clique
    // Need to make sure to not create circular references
    vector<sharedClique> otherChildren;
    otherChildren.insert(otherChildren.begin(), 
                         otherClique->children.begin(),
                         otherClique->children.end());
    for(sharedClique childClique : otherChildren) {
        if(childClique != get_ptr()) {
            childClique->setParent(get_ptr());
        }
    }

    // clique1 parent must now be set to clique2 parent
    setParent(otherClique->parent);

    // After reassigning parent, no node shoud point to this clique
    otherClique->detachParent();
}

void CholeskyEliminationTree::Clique::printClique(ostream& os) {
    os << "Clique: ";
    for(sharedNode node : nodes) {
        os << node->key << " ";
    }
    os << endl;
    os << "   Parent: ";
    if(parent != nullptr) {
        for(sharedNode node : parent->nodes) {
            os << node->key << " ";
        }
    }
    else {
        os << "nullptr";
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
