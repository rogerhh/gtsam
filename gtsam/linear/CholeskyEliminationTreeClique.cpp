/**
* @file    CholeskyEliminationTreeClique.cpp
* @brief   A group of fully connected node in the CholeskyEliminationTree
* @author  Roger Hsiao (rogerhh)
* @date    Mar. 16, 2023
*/

#include "CholeskyEliminationTree.h"
#include <gtsam/linear/CholeskyEliminationTreeClique.h>
#include <gtsam/linear/CholeskyEliminationTreeNode.h>
#include <iostream>
#include <cassert>

using namespace std;

namespace gtsam {

using sharedNode = CholeskyEliminationTree::sharedNode;
using sharedClique = CholeskyEliminationTree::sharedClique;

CholeskyEliminationTree::Clique::Clique(CholeskyEliminationTree* etree_in)
  : etree(etree_in), orderingVersion(etree->orderingVersion_) {}

void CholeskyEliminationTree::Clique::addNode(sharedNode node) {
  nodes.push_back(node);
  node->setClique(shared_from_this());
}

sharedClique CholeskyEliminationTree::Clique::markClique(
    const RemappedKey lowestKey, RemappedKeySet* markedKeys) {

  sharedClique oldParent = parent();
  detachParent();

  if(nodes.size() > 1) {
    // This is not a new clique, check that we have Atb row
    assert(blockIndices.size() >= nodes.size() + 1);
  }

  int i;
  for(i = nodes.size() - 1; i >= 0; i--) {
    markedKeys->insert(nodes[i]->key);
    if(i > 0) {
      // make new clique for node
      sharedClique newClique = make_shared<Clique>(etree);
      newClique->addNode(nodes[i]);
      newClique->marked_ = true;

      assert(gatherSources.size() == 1);
      newClique->gatherSources.push_back(gatherSources.front().split(i));

      if(nodes[i]->key == lowestKey) {
        // The last detached node's clique should have this clique as a child
        setParent(newClique);
        // resize vector to detach all nodes that have new cliques
        nodes.resize(i);
        break;
      }
    }
    else {
      marked_ = true;
      assert(nodes[i]->key == lowestKey);
      assert(nodes[i]->clique() == get_ptr());
      nodes.resize(1);

      // This cannot be a split clique, as clique generated from the previous 
      // if block will not be processed again
      if(!blockIndices.empty()) {
        // If not a new clique
        assert(ownsColumns());

        // This clique now has the same number of columns as the width of the first variable
        // This is taken care of by splitting CliqueColumns

        blockIndices.clear();
      }
    }
  }
  return oldParent;
}

void CholeskyEliminationTree::Clique::reorderClique() {
  orderingVersion = etree->orderingVersion_;

  // Reorder blockIndices
  size_t lowestReorderedIndex = 0;
  vector<RemappedKey> colStructure;
  colStructure.reserve(blockIndices.size());
  for(size_t i = 0; i < blockIndices.size(); i++) {
    colStructure.push_back(get<BLOCK_INDEX_KEY>(blockIndices[i]));
  }
  std::sort(colStructure.begin(), colStructure.end(), etree->orderingLess_);

  bool needReorder = false;
  for(size_t i = 0; i < blockIndices.size() - 1; i++) {
    if(colStructure[i] != get<BLOCK_INDEX_KEY>(blockIndices[i])) {
      needReorder = true;
      break;
    }
  }

  if(needReorder) {
    populateBlockIndices(colStructure);

    assert(ownsColumns());

    // Permute subdiagonal blocks
    if(!gatherSources.empty()) {
      size_t r = height();
      size_t c = width();
      auto matrixSource = std::make_shared<vector<double>>(r * c, 0);
      auto blockIndicesSource = std::make_shared<BlockIndexVector>(blockIndices);

      LocalCliqueColumns newCliqueColumns(matrixSource,
          blockIndicesSource,
          0, cliqueSize());
      newCliqueColumns.addCliqueColumns(gatherSources[0], true);

      gatherSources[0] = newCliqueColumns;
    }
  }
}

sharedClique CholeskyEliminationTree::Clique::parent() {
  sharedClique p = parent_.lock();
  return p;
}

void CholeskyEliminationTree::Clique::findParent() {
  assert(!blockIndices.empty());
  if(blockIndices.size() == 1) {
    // Last row is going to be the root of all cliques
    // This will solve the multiple roots problem
    assert(nodes.front()->key == 0);
    parent_.reset();
    return;
  }

  RemappedKey parentKey = get<BLOCK_INDEX_KEY>(blockIndices.at(cliqueSize()));
  setParent(etree->nodes_[parentKey]->clique());

}

void CholeskyEliminationTree::Clique::reorderAndFindParent() {
  reorderClique();

  findParent();

  assert(parent() != nullptr);
}

void CholeskyEliminationTree::Clique::detachParent() {
  if(parent() != nullptr) {
    assert(parent()->children.find(get_ptr()) != parent()->children.end());
    parent()->children.erase(get_ptr());
    parent_.reset();
  }
}

void CholeskyEliminationTree::Clique::setParent(sharedClique newParent) {
    // Detach from old parent
    detachParent();

    // Set new parent
    parent_ = newParent;
    if(parent() != nullptr) {
        assert(parent()->children.find(get_ptr()) == parent()->children.end());
        parent()->children.insert(get_ptr());
    }
}

shared_ptr<CholeskyEliminationTree::Node> CholeskyEliminationTree::Clique::front() {
    return nodes.front();
}

shared_ptr<CholeskyEliminationTree::Node> CholeskyEliminationTree::Clique::back() {
  return nodes.back();
}

RemappedKey CholeskyEliminationTree::Clique::frontKey() const {
  return nodes.front()->key;
}

RemappedKey CholeskyEliminationTree::Clique::backKey() const {
  return nodes.back()->key;
}

bool CholeskyEliminationTree::Clique::isLastRow() const {
  if(frontKey() == 0) {
    assert(backKey() == 0);
    assert(nodes.size() == 1);
    return true;;
  }
  return false;
}

bool CholeskyEliminationTree::Clique::marked() const {
  return marked_;
}

size_t CholeskyEliminationTree::Clique::cliqueSize() const {
  return nodes.size();
}

size_t CholeskyEliminationTree::Clique::diagonalHeight() const {
  if(isLastRow()) {
    return 1;
  }
  assert(blockIndices.size() > cliqueSize());
  return get<BLOCK_INDEX_ROW>(blockIndices[cliqueSize()]);
}

size_t CholeskyEliminationTree::Clique::subdiagonalHeight() const {
  return get<BLOCK_INDEX_ROW>(blockIndices.back()) - diagonalHeight() + 1;
}

size_t CholeskyEliminationTree::Clique::height() const {
    assert(get<BLOCK_INDEX_KEY>(blockIndices.back()) == 0);
    return get<BLOCK_INDEX_ROW>(blockIndices.back()) + 1;
}

size_t CholeskyEliminationTree::Clique::width() const {
  return diagonalHeight();
}

void CholeskyEliminationTree::Clique::mergeClique(sharedClique otherClique) {
    // Assert the clique is standalone because we generally only merge with lone nodes
    // during symbolic elimination
    assert(parent() == otherClique);
    assert(otherClique->nodes.size() == 1);   
    assert(otherClique->nodes.front()->key != 0);
    assert(otherClique->children.find(get_ptr()) != otherClique->children.end());   

    addNode(otherClique->front());

    // Added clique can only have 1 node, but may have clique children
    // Need to add clique children to merged clique
    // Need to make sure to not create circular references
    // This is for the case of a relaxed supernode. We are disallowing the case for now
    assert(otherClique->children.size() == 1);
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
    assert(otherClique->parent() == nullptr);
    setParent(otherClique->parent());

    // After reassigning parent, no node shoud point to this clique
    otherClique->detachParent();

    // Merge CliqueColumns, which are fragments of the previous Cholesky factor
    mergeGatherSources(otherClique->gatherSources);
}

void CholeskyEliminationTree::Clique::mergeGatherSources(
    const vector<LocalCliqueColumns>& otherGatherSources) {
  // There are 3 cases
  // 1. The 2 cliqueColumns used to belong in the same clique. Merge them
  // 2. The 2 cliqueColumns do not belong in the same clique. Keep them both.
  // 3. One or both of them are new (gatherSources.empty()). Keep whatever is valid.

  if(otherGatherSources.empty()) {
    // Case 3.1. Do nothing
  }
  else if(gatherSources.empty()) {
    // Case 3.2. 
    gatherSources = otherGatherSources;
  }
  else {
    // otherClique should only have a single node, but we can write it to be more flexible
    for(const LocalCliqueColumns& otherColumns : otherGatherSources) {
      // If the other CliqueColumns can be merged into our last CliqueColumn, merge
      // Otherwise just add it to the list of fragments
      if(!gatherSources.back().merge(otherColumns)) {
        gatherSources.push_back(otherColumns);
      }
    }
  }
}

void CholeskyEliminationTree::Clique::mergeColStructure(
    vector<RemappedKey>* parentColStructure) {
  vector<RemappedKey> newColStructure;
  newColStructure.reserve(blockIndices.size() + parentColStructure->size());

  size_t i1 = 0, i2 = cliqueSize();  // merged col structure does not include child diagonal

  while(i1 < parentColStructure->size() || i2 < blockIndices.size()) {
    RemappedKey k1 = 0, k2 = 0; // Default to 0 because the ordering of key 0 is INT_MAX
    if(i1 < parentColStructure->size()) {
      k1 = (*parentColStructure)[i1];
    }
    if(i2 < blockIndices.size()) {
      k2 = get<BLOCK_INDEX_KEY>(blockIndices[i2]);
    }

    if(k1 == k2) {
      newColStructure.push_back(k1);
      i1++;
      i2++;
    }
    else if(etree->orderingLess_(k1, k2)) {
      newColStructure.push_back(k1);
      i1++;
    }
    else {
      assert(etree->orderingLess_(k2, k1));
      newColStructure.push_back(k2);
      i2++;
    }
  }

  *parentColStructure = std::move(newColStructure);
}

void CholeskyEliminationTree::Clique::populateBlockIndices(
    const std::vector<RemappedKey>& colStructure) {
  assert(colStructure.back() == 0);
  blockIndices.clear();

  size_t row = 0;
  for(RemappedKey k : colStructure) {
    size_t height = etree->colWidth(k);
    blockIndices.push_back({k, row, height});
    row += height;
  }
}

bool CholeskyEliminationTree::Clique::hasMarkedAncestor() {
  // only need to check highest key that is not 0 (last row)
  assert(blockIndices.size() > cliqueSize() + 1);
  assert(get<BLOCK_INDEX_KEY>(blockIndices.back()) == 0);

  RemappedKey highestKey = get<BLOCK_INDEX_KEY>(blockIndices[blockIndices.size() - 2]);
  assert(highestKey != 0);

  if(etree->cliques_[highestKey]->marked()) {
    assert(etree->cliques_[highestKey]->status() == RECONSTRUCT);
    return true;
  }
  return false;
}

void CholeskyEliminationTree::Clique::setEditOrReconstruct() {
  assert(0);
}

void CholeskyEliminationTree::Clique::checkEditOrReconstruct(
    MarkedStatus mode, std::vector<RemappedKey>* destCols) {

  destCols->clear();
  destCols->reserve(blockIndices.size());

  for(size_t i = cliqueSize(); i < blockIndices.size() - 1; i++) { 
    // Start from subdiagonal keys and ignore last row
    RemappedKey key = get<BLOCK_INDEX_KEY>(blockIndices[i]);
    assert(key != 0);
    if(etree->nodes_[key]->clique()->status() == mode) {
      // ignore last row
      destCols->push_back(key);
      assert(etree->nodes_[key]->clique()->marked());
    }
  }
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

void CholeskyEliminationTree::Clique::setNodeStatus() {
  // FIXME: delegate this checking to the nodes
  for(sharedNode node : nodes) {
    if(node->status() == NEW && this->status() == EDIT) {
      // Do nothing
    }
    else {
      if(this->status() == EDIT) {
        node->setStatusEdit();
      }
      else if(this->status() == RECONSTRUCT) {
        node->setStatusReconstruct();
      }
      else {
        throw runtime_error("Set invalid node status");
      }
    }
  }
}

void CholeskyEliminationTree::Clique::setNodeStatusMarginalize() {
  for(sharedNode node : nodes) {
    node->setStatusMarginalized();
  }
}

void CholeskyEliminationTree::Clique::setBacksolve(bool backsolve) {
  for(sharedNode node : nodes) {
    node->backsolve = backsolve;
  }
}

bool CholeskyEliminationTree::Clique::needsBacksolve() const {
  assert(get<BLOCK_INDEX_KEY>(blockIndices.back()) == 0);
  // Don't check key 0. Start from the back as higher keys are more likely to 
  // need to backsolve
  for(size_t i = blockIndices.size() - 2; i >= cliqueSize(); i--) {
    RemappedKey key = get<BLOCK_INDEX_KEY>(blockIndices[i]);
    if(etree->nodes_[key]->backsolve) {
      return true;
    }
  }
  return false;
}

sharedClique CholeskyEliminationTree::Clique::splitClique(size_t splitIndex) {
  assert(splitIndex < cliqueSize() && splitIndex > 0);

  sharedClique newClique = make_shared<Clique>(etree);

  for(size_t i = splitIndex; i < cliqueSize(); i++) {
    newClique->addNode(nodes[i]);
  }
  this->nodes.resize(splitIndex);

  newClique->setParent(this->parent());
  this->setParent(newClique);

  assert(this->marked_ == false);
  assert(this->status() == MARGINALIZED);
  assert(this->workspaceIndex == -1);

  size_t startRow = get<BLOCK_INDEX_ROW>(this->blockIndices[splitIndex]);
  for(size_t i = splitIndex; i < blockIndices.size(); i++) {
    newClique->blockIndices.push_back(blockIndices[i]);
    get<BLOCK_INDEX_ROW>(newClique->blockIndices.back()) -= startRow;
  }

  newClique->orderingVersion = this->orderingVersion;

  assert(this->gatherSources.size() == 1);
  assert(this->ownsColumns());

  LocalCliqueColumns newGatherSource = this->gatherSources.front().split(splitIndex);

  newGatherSource.forceOwn();

  newClique->gatherSources.push_back(newGatherSource);

  assert(newClique->ownsColumns());

  for(auto node : this->nodes) {
    assert(node->clique() == this->get_ptr());
  }
  for(auto node : newClique->nodes) {
    assert(node->clique() == newClique);
  }

  return newClique;

}

void CholeskyEliminationTree::Clique::resetAfterCholesky() {
}

void CholeskyEliminationTree::Clique::resetAfterBacksolve() {
  // Reset node member variables
  setBacksolve(false);

  // Reset node member variables
  for(sharedNode node : nodes) {
    node->setStatusUnmarked();
    // node->changedLambdaStructure.clear();
  }

  // Reset member variables
  status_ = UNMARKED;
  marked_ = false;
  workspaceIndex = -1;
}

void CholeskyEliminationTree::Clique::printBlockIndices(ostream& os) {
  for(const auto&[key, row, height] : blockIndices) {
    os << "[" << key << " " << row << " " << height << "] ";
  }
}

void CholeskyEliminationTree::Clique::printClique(ostream& os) {
    os << "Clique: ";
    for(sharedNode node : nodes) {
        os << node->key << " ";
    }
    os << endl;
    os << "   Parent: ";
    if(parent() != nullptr) {
        for(sharedNode node : parent()->nodes) {
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

// Check if we own our columns
bool CholeskyEliminationTree::Clique::ownsColumns() const {

  if(gatherSources.empty()) { return true; }

  if(gatherSources.size() > 1) {
    return false;
  }
  
  return gatherSources.front().ownsData();
}

void CholeskyEliminationTree::Clique::checkInvariant() const {
  for(sharedNode node : nodes) {
    assert(node->clique() == shared_from_this());
  }

  for(int i = 1; i < nodes.size(); i++) {
    assert(etree->orderingLess_(nodes[i - 1]->key, nodes[i]->key));
  }

  for(sharedClique child : children) {
    assert(etree->orderingLess_(child->backKey(), frontKey()));
  }
}

}   // namespace gtsam
