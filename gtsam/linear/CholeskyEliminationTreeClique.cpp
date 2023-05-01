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
    : eTreePtr(eTreePtr_in) {
    gatherSources.push_back(make_tuple(nullptr, 0, 0, 0, 0));
}

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
    if(nodes.size() > 1) {
        // This is not a new clique
        assert(blockIndices.size() >= nodes.size() + 1);
    }

    assert(gatherSources.size() == 1);
    auto&[col_data_ptr, col_data_start, col_row_start, r, c] = gatherSources[0];

    int i;
    for(i = nodes.size() - 1; i >= 0; i--) {
        markedKeys->insert(nodes[i]->key);
        if(i > 0) {
            // make new clique for node
            sharedClique newClique = make_shared<Clique>(eTreePtr);
            newClique->addNode(nodes[i]);
            newClique->marked = true;

            // There is a small problem here, when we split off columns from cliques
            // The new columns shouldn't have the same number of rows as 
            // the top of the column is not used
            size_t newCol = blockIndices[i].second.first;
            newClique->gatherSources.front() = make_tuple(col_data_ptr, 
                                                          newCol * r, 
                                                          newCol, 
                                                          r, blockIndices[i].second.second);
            // newClique->col_data_ptr = col_data_ptr;
            // // New clique will have the same number of rows
            // newClique->col_data_start = newCol * r;
            // newClique->col_row_start = newCol;
            // newClique->r = r;
            // // New clique will have the number of cols as the key
            // newClique->c = blockIndices[i].second.second;

            if(nodes[i]->key == lowestKey) {
                // The last detached node's clique should have this clique as a child
                setParent(newClique);
                // resize vector to detach all nodes that have new cliques
                nodes.resize(i);
                // This clique now has fewer columns than before
                c = blockIndices[i].second.first; // check this
                break;
            }
        }
        else {
            marked = true;
            assert(nodes[i]->key == lowestKey);
            assert(nodes[i]->clique == get_ptr());
            nodes.resize(1);

            // This cannot be a split clique, as clique generated from the previous 
            // if block will not be processed again
            if(!blockIndices.empty()) {
                // If not a new clique
                assert(col_data_ptr != nullptr);
                assert(col_data_start == 0);
                assert(col_row_start == 0);

                // This clique now has the same number of columns as the first variable
                c = blockIndices[0].second.second;
            }
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

size_t CholeskyEliminationTree::Clique::cliqueSize() {
    return nodes.size();
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
    assert(parent == otherClique);
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


    // Merge column matrices. After merging, col_data_start, r, c should be enough information 
    // to gather the existing columns of the clique
    // 3 possibilities
    // 1. Column is just reordered. In that case, this merging does not matter as 
    // we will not use the data. However, we still need to take care of the old memory
    // 2. The parent column is new. In that case, just use the old one
    // 3. The two cliques used to belong to the same clique. In that case, they will have the same col_data_ptr
    //  3.1 However, it is also possible that there is a reordering but the two nodes are still in the same clique. In this case, the old col will still be downsized the same amount
    // I.e. col_data_ptr should be the same for all 3 cases. r and c need to be merged in case 3

    assert(otherClique->gatherSources.size() == 1);

    auto&[col_data_ptr1, col_data_start1, col_row_start1, r1, c1] = gatherSources.back();
    auto&[col_data_ptr2, col_data_start2, col_row_start2, r2, c2] = otherClique->gatherSources.back();
    
    if(col_data_ptr1 != nullptr && col_data_ptr1 == col_data_ptr2) {
        // This is case 3
        assert(r1 == r2);
        c1 += c2;
    }
    else {
        // If col_data_ptr == nullptr
        //   In this case, either both are new or there is a reordering. Might need to deallocate parent
        // If col_data_ptr != otherClique->col_data_ptr
        //   There is a reordering, might need to deallocate parent
        //   (New scheme: don't deallocate old but keep it around)
        if(col_data_ptr2) {
            gatherSources.push_back(make_tuple(col_data_ptr2, col_data_start2, 
                                               col_row_start2, r2, c2));
            /*
            size_t otherSize = col_data_ptr2->size();
            size_t otherR = r;
            size_t otherC = otherClique->c;
            assert(otherSize >= r2 * c2);
            col_data_ptr2->resize(otherSize - r2 * otherC);
            */

        }
        else {
        }
    
    }
}

bool CholeskyEliminationTree::Clique::reorderColumn(BlockIndexVector& newBlockIndices) {
    assert(newBlockIndices.size() == blockIndices.size());

    size_t i = 0;
    // The nodes in the clique should remain the same
    for(i = 0; i < cliqueSize(); i++) {
        assert(newBlockIndices[i] == blockIndices[i]);
    }

    size_t lowestReorderedIndex = -1;
    for(; i < blockIndices.size(); i++) {
        if(newBlockIndices[i].first != blockIndices[i].first) {
            lowestReorderedIndex = i;
            break;
        }
    }

    if(i == blockIndices.size()) {
        return false;
    }

    assert(blockIndices[lowestReorderedIndex].second.first 
            == newBlockIndices[lowestReorderedIndex].second.first);
    size_t firstRow = blockIndices[lowestReorderedIndex].second.first;

    unordered_map<Key, size_t> keyRowMap;

    for(size_t i = lowestReorderedIndex; i < blockIndices.size(); i++) {
        Key key = blockIndices[i].first;
        size_t row = blockIndices[i].second.first;
        keyRowMap.insert({key, row});
    }

    assert(gatherSources.size() == 1);
    auto&[col_data_ptr, col_data_start, col_row_start, r, c] = gatherSources[0];

    // if(firstRow * 2 < r) {
    if(false) {
        // Update partially
        assert(0);
    
    } 
    else {
        // update all matrix
        shared_ptr<vector<double>> new_col_data_ptr = make_shared<vector<double>>(r * c);
        Eigen::Map<ColMajorMatrix> new_m(new_col_data_ptr->data(), r, c);
        Eigen::Map<ColMajorMatrix> old_m(col_data_ptr->data(), r, c);

        Eigen::Block<Eigen::Map<ColMajorMatrix>>(new_m, 0, 0, firstRow, c)
            = Eigen::Block<Eigen::Map<ColMajorMatrix>>(old_m, 0, 0, firstRow, c);

        for(size_t i = lowestReorderedIndex; i < blockIndices.size(); i++) {
            Key newKey = newBlockIndices[i].first;
            size_t newRow = newBlockIndices[i].second.first;
            size_t width = newBlockIndices[i].second.second;
            size_t oldRow = keyRowMap.at(newKey);

            Eigen::Block<Eigen::Map<ColMajorMatrix>>(new_m, newRow, 0, width, c)
                = Eigen::Block<Eigen::Map<ColMajorMatrix>>(old_m, oldRow, 0, width, c);
        }

        col_data_ptr = new_col_data_ptr;

    }

    blockIndices = std::move(newBlockIndices);
    
    return true;
}

// void CholeskyEliminationTree::Clique::deallocateExcessCols() {
//     assert(!marked);
//     assert(temp_data_start == 0);
// 
//     if(r * c != col_data.size()) {
//         assert(parent != nullptr);
//         assert(parent->temp_data == &col_data);
//         parent->temp_data = nullptr;
//         col_data.resize(r * c);
//     }
//     else {
//         assert(parent != nullptr);
//         assert(parent->temp_data != &col_data);
//     }
// }

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

std::ostream& operator<<(std::ostream& os, const CholeskyEliminationTree::Clique& clique) {
    for(size_t i = 0; i < clique.nodes.size(); i++) {
        os << clique.nodes[i]->key;
        if(i != clique.nodes.size() - 1) {
            os << " ";
        }
    }
    return os;
}


}   // namespace gtsam
