/**
* @file    CholeskyEliminationTreeClique.cpp
* @brief   A group of fully connected node in the CholeskyEliminationTree
* @author  Roger Hsiao (rogerhh)
* @date    Mar. 16, 2023
*/

#include "CholeskyEliminationTree.h"
#include "gtsam/linear/CholeskyEliminationTreeTypes.h"
#include <gtsam/linear/CholeskyEliminationTreeClique.h>
#include <gtsam/linear/CholeskyEliminationTreeFactorWrapper.h>
#include <gtsam/linear/CholeskyEliminationTreeNode.h>
#include <iostream>
#include <cassert>
#include <chrono>

#include <gtsam/linear/model.h>

using namespace std;

namespace gtsam {

using sharedNode = CholeskyEliminationTree::sharedNode;
using sharedClique = CholeskyEliminationTree::sharedClique;

uint64_t CholeskyEliminationTree::Clique::merge1_cycles = 0;
uint64_t CholeskyEliminationTree::Clique::merge2_cycles = 0;
uint64_t CholeskyEliminationTree::Clique::merge3_cycles = 0;
uint64_t CholeskyEliminationTree::Clique::merge4_cycles = 0;
uint64_t CholeskyEliminationTree::Clique::merge5_cycles = 0;
uint64_t CholeskyEliminationTree::Clique::merge6_cycles = 0;
uint64_t CholeskyEliminationTree::Clique::gather1_cycles = 0;
uint64_t CholeskyEliminationTree::Clique::gather2_cycles = 0;
uint64_t CholeskyEliminationTree::Clique::gather3_cycles = 0;
uint64_t CholeskyEliminationTree::Clique::model1_cycles = 0;
uint64_t CholeskyEliminationTree::Clique::model2_cycles = 0;

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
      auto matrixSource = std::make_shared<vector<GEMMINI_TYPE>>(r * c, 0);
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

void CholeskyEliminationTree::Clique::setParent2(sharedClique newParent) {
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

void CholeskyEliminationTree::Clique::mergeClique2(sharedClique childClique) {
    // Merge children into ourself
    assert(childClique->parent() == get_ptr());
    assert(nodes.front()->key != 0);
    assert(children.find(childClique) != children.end());

    auto merge1_start = chrono::high_resolution_clock::now();

    // First merge blockIndices
    int i1 = 0, i2 = 0;
    BlockIndexVector thisBlockIndices = std::move(this->blockIndices);
    const BlockIndexVector& childBlockIndices = childClique->blockIndices;
    this->blockIndices.reserve(thisBlockIndices.size() + childBlockIndices.size());
    int curRow = 0;
    auto merge1_end = chrono::high_resolution_clock::now();

    while(i1 < childBlockIndices.size() && i2 < thisBlockIndices.size()) {
      auto&[k1, _1, h1] = childBlockIndices[i1];
      auto&[k2, _2, h2] = thisBlockIndices[i2];
      if(k1 == k2) {
        this->blockIndices.push_back({k1, curRow, h1});
        curRow += h1;
        i1++;
        i2++;
      }
      else if(etree->orderingLess_(k1, k2)) {
        this->blockIndices.push_back({k1, curRow, h1});
        curRow += h1;
        i1++;
      }
      else {
        this->blockIndices.push_back({k2, curRow, h2});
        curRow += h2;
        i2++;
      }
    }

    while(i1 < childBlockIndices.size()) {
      auto&[k1, _1, h1] = childBlockIndices[i1];
      this->blockIndices.push_back({k1, curRow, h1});
      curRow += h1;
      i1++;
    }

    while(i2 < thisBlockIndices.size()) {
      auto&[k2, _2, h2] = thisBlockIndices[i2];
      this->blockIndices.push_back({k2, curRow, h2});
      curRow += h2;
      i2++;
    }

    auto merge2_end = chrono::high_resolution_clock::now();

    int cliqueSize = childClique->cliqueSize() + this->cliqueSize();

    this->nodes.resize(cliqueSize);
    for(int i = 0; i < cliqueSize; i++) {
      RemappedKey k = get<BLOCK_INDEX_KEY>(this->blockIndices[i]);
      this->nodes[i] = etree->nodes_[k];
      etree->cliques_[k] = get_ptr();
    }
    auto merge3_end = chrono::high_resolution_clock::now();

    // Detach child cliques parent
    // At this point, nothing should point to child anymore
    childClique->detachParent();

    auto merge4_end = chrono::high_resolution_clock::now();

    // Added clique may have children. We need to adopt the child cliques children
    for(sharedClique clique : childClique->children) {
        clique->setParent2(get_ptr());
    }
    // vector<sharedClique> childChildren;
    // childChildren.insert(childChildren.begin(), 
    //                      childClique->children.begin(),
    //                      childClique->children.end());
    // for(sharedClique childClique : childChildren) {
    //     childClique->setParent(get_ptr());
    // }

    auto merge5_end = chrono::high_resolution_clock::now();

    // Merge CliqueColumns, which are fragments of the previous Cholesky factor
    mergeGatherSources2(childClique->gatherSources);

    auto merge6_end = chrono::high_resolution_clock::now();

    merge1_cycles += chrono::duration_cast<chrono::nanoseconds>(merge1_end - merge1_start).count();
    merge2_cycles += chrono::duration_cast<chrono::nanoseconds>(merge2_end - merge1_end).count();
    merge3_cycles += chrono::duration_cast<chrono::nanoseconds>(merge3_end - merge2_end).count();
    merge4_cycles += chrono::duration_cast<chrono::nanoseconds>(merge4_end - merge3_end).count();
    merge5_cycles += chrono::duration_cast<chrono::nanoseconds>(merge5_end - merge4_end).count();
    merge6_cycles += chrono::duration_cast<chrono::nanoseconds>(merge6_end - merge5_end).count();

}

void CholeskyEliminationTree::Clique::mergeGatherSources2(
    vector<LocalCliqueColumns>& childGatherSources) {
  // There are 3 cases
  // 1. The 2 cliqueColumns used to belong in the same clique. Merge them
  // 2. The 2 cliqueColumns do not belong in the same clique. Keep them both.
  // 3. One or both of them are new (gatherSources.empty()). Keep whatever is valid.

  if(childGatherSources.empty()) {
    // Case 3.1. Do nothing
  }
  else if(this->gatherSources.empty()) {
    // Case 3.2. 
    auto gather3_start = chrono::high_resolution_clock::now();
    this->gatherSources = std::move(childGatherSources);
    auto gather3_end = chrono::high_resolution_clock::now();
    gather3_cycles += chrono::duration_cast<chrono::nanoseconds>(gather3_end - gather3_start).count();
  }
  else {
    // otherClique should only have a single node, but we can write it to be more flexible
    auto gather1_start = chrono::high_resolution_clock::now();
    std::swap(childGatherSources, this->gatherSources);
    auto gather1_end = chrono::high_resolution_clock::now();
    for(const LocalCliqueColumns& columns : childGatherSources) {
      // If our CliqueColumns can be merged into the child's last CliqueColumn, merge
      // Otherwise just add it to the list of fragments
      if(!this->gatherSources.back().merge(columns)) {
        this->gatherSources.push_back(columns);
      }
    }
    auto gather2_end = chrono::high_resolution_clock::now();
    gather1_cycles += chrono::duration_cast<chrono::nanoseconds>(gather1_end - gather1_start).count();
    gather2_cycles += chrono::duration_cast<chrono::nanoseconds>(gather2_end - gather1_end).count();
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

void CholeskyEliminationTree::Clique::extendBlockIndices(
    const std::vector<RemappedKey>& colStructure) {
  const size_t cliqueSize = this->cliqueSize();

  assert(colStructure.back() == 0);
  assert(get<BLOCK_INDEX_KEY>(blockIndices[cliqueSize - 1]) == colStructure[0]);

  size_t row = get<BLOCK_INDEX_ROW>(blockIndices[cliqueSize]);
  blockIndices.resize(cliqueSize - 1 + colStructure.size());

  for(size_t i = 1; i < colStructure.size(); i++) {
    RemappedKey k = colStructure[i];
    size_t height = etree->colWidth(k);
    blockIndices[cliqueSize - 1 + i] = {k, row, height};
    row += height;
  }

}

bool CholeskyEliminationTree::Clique::hasMarkedAncestor() {
  // only need to check highest key that is not 0 (last row)
  assert(blockIndices.size() >= cliqueSize() + 1);
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

void CholeskyEliminationTree::Clique::deleteClique() {
  this->detachParent();

  if(!children.empty()) {
    cout << "Deleted clique should not have children!" << endl;
    exit(1);
  }
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


void CholeskyEliminationTree::Clique::resetCost() {
    costStatus = COST_UNMARKED;
    nextCostStatus = COST_UNMARKED;
    markedCost = -1;
    fixedCost = -1;
    backsolveCost = -1;
}

const double cost_mplier = 1.25;

int64_t CholeskyEliminationTree::Clique::computeCostMarked(int num_threads) {

  if(markedCost >= 0) { return markedCost * cost_mplier; }


  auto model1_start = chrono::high_resolution_clock::now();
  int maxFactorHeight = 0;
  int maxFactorWidth = 0;
  int num_factors = 0;
  for(sharedNode node : nodes) {
    for(sharedFactorWrapper factorWrapper : node->factors) {
      int r, c;

      if(node->key != factorWrapper->lowestKey()) { continue; }

      if(factorWrapper->status() == LINEARIZED || factorWrapper->status() == LINEAR) {
        r = factorWrapper->getCachedMatrix().rows();
        c = factorWrapper->getCachedMatrix().cols();
      }
      else {
        r = factorWrapper->nonlinearFactor()->dim();
        c = get<BLOCK_INDEX_ROW>(factorWrapper->blockIndices().back()) + get<BLOCK_INDEX_HEIGHT>(factorWrapper->blockIndices().back());
      }

      if(r > maxFactorHeight) { maxFactorHeight = r; }
      if(c > maxFactorWidth) { maxFactorWidth = c; }

      num_factors++;
    }
  }
  auto model1_end = chrono::high_resolution_clock::now();

  int cliqueWidth = this->width();
  int cliqueHeight = this->height();
  int chol_tile_len = 8;
  int gemm_tile_len = 32;
  bool next_memset;
  int next_node_height;

  if(this->parent() == nullptr || this->parent() == etree->root_) {
    next_memset = false;
    next_node_height = 0;
  }
  else {
    next_memset = true;
    next_node_height = this->parent()->height();
  }

  // Factor is transposed
  int64_t pred_restore = predict_restore(cliqueWidth, cliqueHeight);
  int64_t pred_AtA = predict_AtA_reorder(maxFactorWidth, maxFactorHeight, num_factors);
  int64_t pred_chol = predict_cholesky(cliqueWidth, cliqueHeight, chol_tile_len, gemm_tile_len, next_memset, next_node_height);
  int64_t pred_add = predict_node_add(cliqueHeight, cliqueWidth);

  int64_t restore_overhead = 0;
  int64_t AtA_overhead = 2500;
  int64_t chol_overhead = 0;
  int64_t add_overhead = 1000;

  auto model2_end = chrono::high_resolution_clock::now();

  // This is used for single node profiling
  restoreCost = pred_restore + restore_overhead;
  symCost = this->nodes.size() * SYM_COST_REORDER;
  AtACost = AtA_overhead + pred_AtA;
  cholCost = chol_overhead + pred_chol;
  addCost = add_overhead + pred_add;

  if(!this->parallelizable) {
    markedCost = restoreCost + symCost + AtACost + cholCost + addCost;
  }
  else {
    markedCost = (pred_restore + pred_AtA + pred_chol) / num_threads 
                + (symCost + pred_add + AtA_overhead + chol_overhead + add_overhead);
  }

  // cout << "Clique " << *this << " ata cost = " << markedCost << " markedCost = " << fixedCost << endl;
  
  model1_cycles += chrono::duration_cast<chrono::nanoseconds>(model1_end - model1_start).count();
  model2_cycles += chrono::duration_cast<chrono::nanoseconds>(model2_end - model1_end).count();

  return markedCost * cost_mplier;

}


int64_t CholeskyEliminationTree::Clique::computeCostFixed(int num_threads) {

  int maxFactorHeight = 0;
  int maxFactorWidth = 0;
  int num_factors = 0;

  for(sharedNode node : nodes) {
    for(sharedFactorWrapper factorWrapper : node->factors) {

      if(node->key != factorWrapper->lowestKey()) { continue; }

      // Only consider factor if a key in it is marked
      bool flag = false;
      for(RemappedKey remappedKey : factorWrapper->remappedKeys()) {
        if(etree->cliques_[remappedKey]->nextCostStatus == COST_MARKED) {
          flag = true;
          break;
        }
      }

      // if(!flag) { continue; }

      int r, c;

      if(factorWrapper->status() == LINEARIZED || factorWrapper->status() == LINEAR) {
        r = factorWrapper->getCachedMatrix().rows();
        c = factorWrapper->getCachedMatrix().cols();
      }
      else {
        r = factorWrapper->nonlinearFactor()->dim();
        c = get<BLOCK_INDEX_ROW>(factorWrapper->blockIndices().back()) + get<BLOCK_INDEX_HEIGHT>(factorWrapper->blockIndices().back());
      }

      if(r > maxFactorHeight) { maxFactorHeight = r; }
      if(c > maxFactorWidth) { maxFactorWidth = c; }

      num_factors++;
    }
  }

  int cliqueWidth = this->width();
  int cliqueHeight = this->height();
  int chol_tile_len = 8;
  int gemm_tile_len = 32;
  bool next_memset;
  int next_node_height;

  if(this->parent() == nullptr || this->parent() == etree->root_) {
    next_memset = false;
    next_node_height = 0;
  }
  else {
    next_memset = true;
    next_node_height = this->parent()->height();
  }

  // Factor is transposed
  int64_t pred_AtA = predict_AtA_reorder(maxFactorWidth, maxFactorHeight, num_factors);
  int64_t pred_syrk = predict_syrk(cliqueWidth, cliqueHeight, next_memset, next_node_height);
  int64_t pred_add = predict_node_add(cliqueHeight, cliqueWidth);

  int64_t AtA_min = 2000;
  int64_t AtA_overhead = pred_AtA < AtA_min? AtA_min : 0;
  int64_t syrk_min = 2000;
  int64_t chol_overhead = pred_syrk < syrk_min? syrk_min : 0;
  int64_t add_overhead = 1000;

  // This is used for single node profiling
  symCost = this->nodes.size() * SYM_COST_REORDER;
  AtACost = AtA_overhead + pred_AtA;
  syrkCost = chol_overhead + pred_syrk;
  addCost = add_overhead + pred_add;

  if(!this->parallelizable) {
    fixedCost = AtACost + syrkCost + addCost + symCost;
  }
  else {
    fixedCost = (pred_AtA + pred_syrk) / num_threads + (symCost + pred_add + AtA_overhead + chol_overhead + add_overhead);
  }

  // cout << "Clique " << *this << " ata cost = " << AtACost << " fixedCost = " << fixedCost << endl;

  return fixedCost * cost_mplier;

}

int64_t CholeskyEliminationTree::Clique::computeCostBacksolve(int num_threads) {

  int cliqueWidth = this->width();
  int cliqueHeight = this->height();
  int backsolve_block_len = 4;

  int64_t backsolve_overhead = 0;
  int64_t pred_backsolve = predict_backsolve(cliqueWidth, cliqueHeight, backsolve_block_len);

  if(!this->parallelizable) {
    backsolveCost = backsolve_overhead + pred_backsolve;
  }
  else {
    backsolveCost = (pred_backsolve) / num_threads + backsolve_overhead;
  }

  return backsolveCost * cost_mplier;

}

}   // namespace gtsam
