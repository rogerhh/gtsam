/**
* @file    CholeskyEliminationTree.h
* @brief   Elimination tree structure to perform symbolic factorization and Cholesky factorization
* @author  Roger Hsiao (rogerhh)
* @date    Feb. 8, 2023
*/

#include "gtsam/nonlinear/NonlinearFactor.h"
#include <gtsam/base/types.h>
#include <gtsam/inference/Ordering.h>
#include <gtsam/linear/CholeskyEliminationTree.h>
#include <gtsam/linear/CholeskyEliminationTreeNode.h>
#include <gtsam/linear/CholeskyEliminationTreeClique.h>
#include <gtsam/linear/CholeskyEliminationTreeFactorWrapper.h>
#include <gtsam/linear/CholeskyEliminationTreeTypes.h>
#include <gtsam/linear/CliqueColumns.h>
#include <gtsam/3rdparty/CCOLAMD/Include/ccolamd.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <stdexcept>

using namespace std;

namespace gtsam {

// Convenience function
inline Eigen::Block<Eigen::Map<ColMajorMatrix>> block(
    Eigen::Map<ColMajorMatrix> m, 
    size_t r, size_t c, size_t h, size_t w) {
    return Eigen::Block<Eigen::Map<ColMajorMatrix>>(m, r, c, h, w);
}

CholeskyEliminationTree::CholeskyEliminationTree() : orderingLess_(this) {
  // Add new node for the last row. Technically should be connected to everything
  RemappedKey key = addRemapKey(-1);
  assert(key == 0);
  addNewNode(key, 1); 
}

void CholeskyEliminationTree::addVariables(const Values& newTheta) {
  // cout << "[CholeskyEliminationTree] addVariables() " << newTheta.size() << endl;
  for(const auto& keyValPair : newTheta) {
    const Key& unmappedKey = keyValPair.key;
    const Value& val = keyValPair.value;
    const size_t dim = val.dim();

    // First remap the key to start from 1
    RemappedKey key = addRemapKey(unmappedKey);
    addNewNode(key, dim);
    
  }
} 

// Mark directly changed keys and keys that we explicitly want to update (extraKeys)
// observedKeys are the keys associated with the new factors
void CholeskyEliminationTree::markAffectedKeys(
  const NonlinearFactorGraph& nonlinearFactors,
  const FactorIndices& newFactorIndices,
  const KeySet& relinKeys, 
  const ISAM2UpdateParams& updateParams,
  RemappedKeySet* affectedKeys,
  RemappedKeySet* observedKeys) {

  // cout << "[CholeskyEliminationTree] markAffectedKeys()" << endl;
  affectedKeys->clear();
  observedKeys->clear();

  // RelinKeys should be processed before we add in factors because we only need to
  // relinearize old factors
  for(const Key unmappedRelinKey : relinKeys) {
    RemappedKey relinKey = getRemapKey(unmappedRelinKey);
    sharedNode relinNode = nodes_[relinKey];
    assert(relinNode->status() != NEW);

    // Ignore all relin nodes that are MARGINALIZED (fixed)
    if(relinNode->status() == MARGINALIZED) {
      continue;
    }

    for(sharedFactorWrapper factorWrapper : relinNode->factors) {
      if(factorWrapper->status() == LINEARIZED) {
        factorWrapper->setStatusRelinearize();
        factorWrapper->markAffectedKeys(affectedKeys);
      }
    }
  }
  
  for(const FactorIndex newFactorIndex : newFactorIndices) {
    // newFactorIndex starts from 0, and does not necessarily correspond to total factor index
    BlockIndexVector blockIndices;
    sharedFactor factor = nonlinearFactors[newFactorIndex];
    size_t factorIndex = factors_.size();
    sharedFactorWrapper factorWrapper = std::make_shared<FactorWrapper>(
                                          factorIndex, factor, nullptr, this);

    // Remap old factorIndex to real factor index
    // Removing factors uses old factorIndex
    if(newFactorIndex >= factorIndexTransformMap_.size()) {
      for(size_t i = factorIndexTransformMap_.size(); i <= newFactorIndex; i++) {
        factorIndexTransformMap_.push_back(-1);
      }
    }
    factorIndexTransformMap_[newFactorIndex] = factorIndex;

    // If factor involves variables that are marginalized, ignore factor
    if(factorWrapper->hasMarginalizedKeys()) {
      assert(0);
    }

    factors_.push_back(factorWrapper);

    for(RemappedKey k : factorWrapper->remappedKeys()) {
      // k is in every factor, we shouldn't have to store it
      if(k != 0) {
        sharedNode node = nodes_[k];
        node->addFactor(factorWrapper);
      }
    }

    factorWrapper->markAffectedKeys(affectedKeys);
    factorWrapper->markAffectedKeys(observedKeys);
  }
  
  for(const FactorIndex removeFactorIndex : updateParams.removeFactorIndices) {
    // Removing factors remove the old factor indices, disregarding all the marginal
    // and linear factors we add in
    size_t realFactorIndex = factorIndexTransformMap_[removeFactorIndex];
    sharedFactorWrapper factorWrapper = factors_.at(realFactorIndex);
    
    // If factor involves variables that are marginalized, it cannot be removed
    if(factorWrapper->hasMarginalizedKeys()) {
      assert(0);
    }

    factorWrapper->setStatusRemoving();
    factorWrapper->markAffectedKeys(affectedKeys);
  }

  // We don't have to worry about noRelinKeys for now as it will be removed from relinKeys
  
  // We have to manually mark extraReelimKeys as affectedKeys though
  if(updateParams.extraReelimKeys) {
    for(Key unmappedExtraReelimKey : updateParams.extraReelimKeys.get()) {
      RemappedKey extraReelimKey = getRemapKey(unmappedExtraReelimKey);
      affectedKeys->insert(extraReelimKey);
    }
  }

  if(updateParams.newAffectedKeys) {
    for(auto&[factorIndex, newUnmappedKeys] : updateParams.newAffectedKeys.get()) {
      vector<RemappedKey> newKeys;
      for(Key newUnmappedKey : newUnmappedKeys) {
        newKeys.push_back(getRemapKey(newUnmappedKey));
      }
      sharedFactorWrapper factorWrapper = factors_[factorIndex];
      factors_[factorIndex]->addNewKeys(newKeys);
      factors_[factorIndex]->markAffectedKeys(affectedKeys);
      factors_[factorIndex]->markAffectedKeys(observedKeys);
    }
  }

  // DEBUG
  for(RemappedKey key : *affectedKeys) {
    if(!(nodes_[key]->status() == UNMARKED || nodes_[key]->status() == NEW)) {
      cout << "Node " << key << " status = " << nodes_[key]->status() << endl;
    }
    assert(nodes_[key]->status() == UNMARKED || nodes_[key]->status() == NEW);
  }
  // DEBUG END

  // // DEBUG
  // cout << "Observed keys: ";
  // for(RemappedKey key : *observedKeys) {
  //   cout << key << " ";
  // }
  // cout << endl;

  // cout << "RelinKeys: ";
  // for(const Key relinKey : relinKeys) {
  //     cout << getRemapKey(relinKey) << " ";
  // }
  // cout << endl;

  // cout << "Affected keys: ";
  // for(RemappedKey key : *affectedKeys) {
  //   cout << key << " ";
  // }
  // cout << endl;
  //  // DEBUG END

}

void CholeskyEliminationTree::markAncestors(
    const RemappedKeySet& affectedKeys, 
    RemappedKeySet* markedKeys) {
  // cout << "[CholeskyEliminationTree] markAncestors()" << endl;
  for(const RemappedKey key : affectedKeys) {
      markKey(key, markedKeys);
  }
  for(Key k : *markedKeys) {  // marked keys are always backsolved
    nodes_[k]->backsolve = true;
  }
}

void CholeskyEliminationTree::markKey(const RemappedKey key, RemappedKeySet* markedKeys) {
  // Not sure why this is needed
  if(nodes_[key]->clique()->marked()) {
    // Node is already processed
    return;
  }

  // cout << "[CholeskyEliminationTree] markKey " << key << " in clique " << *(nodes_[key]->clique()) << endl;
  sharedNode node = nodes_[key];

  sharedClique curClique = node->clique();
  RemappedKey curKey = key;
  do {
    // if(curClique->status() != UNMARKED) {
    //   cout << "mark clique " << *curClique << " is not UNMARKED. status = " << curClique->status() << endl;
    // }
    assert(curClique->status() == UNMARKED);

    if(curClique->orderingVersion != orderingVersion_) {
      curClique->reorderClique();
    }

    curClique = curClique->markClique(curKey, markedKeys);
    if(curClique) {
      curKey = curClique->frontKey();
    }
  } while(curClique != nullptr);
}

void CholeskyEliminationTree::symbolicElimination(const RemappedKeySet& markedKeys) {
  // cout << "[CholeskyEliminationTree] symbolicElimination()" << endl;

  root_ = nullptr;

  // cout << "markedKeys: ";
  // for(RemappedKey key : markedKeys) {
  //   cout << key << " ";
  // }
  // cout << endl;

  // Sort marked keys based on variable ordering
  vector<RemappedKey> sortedMarkedKeys;
  sortedMarkedKeys.reserve(markedKeys.size());
  sortedMarkedKeys.insert(sortedMarkedKeys.begin(), markedKeys.begin(), markedKeys.end());
  sort(sortedMarkedKeys.begin(), sortedMarkedKeys.end(), orderingLess_);

  for(const RemappedKey key : sortedMarkedKeys) {
    symbolicEliminateKey(key);
  }
  assert(root_);

  if(postOrder_) {
    // sortedMarkedKeys in this case will be exactly the same as partialOrdering
    postReordering(sortedMarkedKeys);
  }

  allocateStack();

#ifdef DEBUG
  checkInvariant_afterSymbolic();
#endif
}

void CholeskyEliminationTree::symbolicEliminateKey(const RemappedKey key) {
  // cout << "[CholeskyEliminationTree] symbolicEliminateKey: " << key << endl;

  sharedNode node = nodes_[key];
  sharedClique clique = node->clique();

  assert(clique->orderingVersion == orderingVersion_);
  assert(clique->nodes.size() == 1);
  assert(clique->marked() == true);

  // Add keys induced by raw factors but only keys that are higher than this key
  vector<RemappedKey> colStructure;
  for(const auto& [otherKey, count] : node->lambdaStructure) {
    if(!orderingLess_(otherKey, node->key)) {
      colStructure.push_back(otherKey);
    }
  }
  std::sort(colStructure.begin(), colStructure.end(), orderingLess_);
  // colStructure.insert(colStructure.end(),
  //                     node->lambdaStructure.find(node->key),
  //                     node->lambdaStructure.end());

  // Merge column structure from children cliques
  for(sharedClique childClique : clique->children) {
      if(childClique->orderingVersion != orderingVersion_) {
          childClique->reorderClique();
      }

      childClique->mergeColStructure(&colStructure);
  }

  // cout << "col structure: ";
  // for(Key k : colStructure) {
  //     cout << k << " ";
  // }
  // cout << endl;

  bool mergeFlag = false;
  if(clique->nodes.front()->key != 0 &&
      clique->children.size() == 1) {
    // Clique 0 cannot be merged and will always be root
    sharedClique childClique = *(clique->children.begin());
    if(childClique->marked() 
        && colStructure.size() == childClique->blockIndices.size() - childClique->cliqueSize()) {
      childClique->mergeClique(clique);
      // Need to update clique pointer to current clique, which is the child clique
      assert(clique.unique());
      clique.reset();
      clique = childClique;
      mergeFlag = true;
    }
  }

  // Need to popoluate our own blockIndices
  if(!mergeFlag) {
    clique->populateBlockIndices(colStructure);
  }

  // Find parent after merging cliques
  clique->findParent();

  // Set root after merging
  if(clique->parent() == nullptr) {
    assert(clique->nodes.size() == 1);
    assert(clique->nodes.front()->key == 0);

    if(root_ != nullptr) {
      // We can actually solve the mltiple root issue by recognizing 
      // the node 0 (Atb row) as a node
      throw runtime_error("Multiple roots!");
    }
    root_ = clique;
  }
}

void CholeskyEliminationTree::remapConstrainedKeys(
    const boost::optional<FastMap<Key, int>>& unmappedConstrainedKeys,
    std::map<RemappedKey, int>* constrainedKeys) {
  constrainedKeys->clear();
  if(unmappedConstrainedKeys) {
    for(const auto&[unmappedKey, group] : *unmappedConstrainedKeys) {
      constrainedKeys->insert({getRemapKey(unmappedKey), group});
    }
  }
}

void CholeskyEliminationTree::constructCSCMatrix(
    const unordered_map<RemappedKey, int>& markedKeysIndex,
    const vector<RemappedKey>& markedKeysVector,
    int* nEntries,
    int* nVars,
    int* nFactors,
    vector<int>* A,
    vector<int>* p) {

  assert(markedKeysVector.back() != 0);

  vector<vector<int>> csc_matrix(markedKeysVector.size());
  csc_matrix = vector<vector<int>>(markedKeysIndex.size());

  unordered_map<sharedFactorWrapper, size_t> rawFactors;

  size_t entriesCount = 0;

  for(const RemappedKey& key : markedKeysVector) {
    sharedNode node = nodes_[key];
    for(sharedFactorWrapper factor : node->factors) {
      if(factor->status() != REMOVING || factor->status() != REMOVED) {
        auto it = rawFactors.insert({factor, rawFactors.size()}).first;
        size_t factorIndex = it->second;
        csc_matrix[markedKeysIndex.at(key)].push_back(factorIndex);
        entriesCount++;
      }
    }
  }
  size_t factorCount = rawFactors.size();

  // Deal with induced factors, each child contributes an induced factor
  for(const RemappedKey& key : markedKeysVector) {
    sharedClique clique = cliques_[key];
    for(sharedClique child : clique->children) {
      assert(!child->marked());
      for(size_t i = child->cliqueSize(); i < child->blockIndices.size() - 1; i++) {
        RemappedKey key = get<BLOCK_INDEX_KEY>(child->blockIndices[i]);
        auto it = markedKeysIndex.find(key);
        assert(it != markedKeysIndex.end());
        csc_matrix[it->second].push_back(factorCount);
        entriesCount++;
      }
      factorCount++;
    }
  }

  p->clear();
  p->reserve(markedKeysVector.size());
  p->push_back(0);

  A->clear();
  A->reserve(entriesCount);

  int count = 0;
  for(const vector<int>& A_col : csc_matrix) {
    // TODO: Check if factors are sorted
    for(int row : A_col) {
      A->push_back(row);
      count++;
    }
    p->push_back(count);
  }

  *nVars = markedKeysVector.size();
  *nFactors = factorCount;
  *nEntries = count;
}

void CholeskyEliminationTree::getPartialReordering(
    const RemappedKeySet& markedKeys, 
    const RemappedKeySet& observedKeys,
    const ISAM2UpdateParams& params,
    vector<RemappedKey>* partialOrdering) {

  // cout << "[CholeskyEliminationTree] getPartialReordering()" << endl;

  unordered_map<RemappedKey, int> markedKeysIndex;
  vector<RemappedKey> markedKeysVector;
  for(RemappedKey key : markedKeys) {
    markedKeysVector.push_back(key);
  }

  std::sort(markedKeysVector.begin(), markedKeysVector.end(), orderingLess_);
  assert(markedKeysVector.back() == 0);

  // No need to add key 0
  markedKeysVector.pop_back();

  for(RemappedKey key : markedKeysVector) {
    markedKeysIndex.insert({key, markedKeysIndex.size()});
  }

  int nEntries, nVars, nFactors;
  vector<int> A, p;
  constructCSCMatrix(markedKeysIndex, markedKeysVector, &nEntries, &nVars, &nFactors, &A, &p);
  const size_t Alen = ccolamd_recommended(nEntries, nFactors, nVars);
  A.resize(Alen);

  // Construct constraints
  map<int, vector<RemappedKey>> constrainedGroupBuckets;
  for(RemappedKey key : markedKeys) {
    bool found = false;
    // If group is already defined by params.constrainedKeys
    if(params.constrainedKeys) {
      Key unmappedKey = unmapKey(key);
      auto it = params.constrainedKeys->find(unmappedKey);
      if(it != params.constrainedKeys->end()) {
        constrainedGroupBuckets[it->second].push_back(key);
        found = true;
      }
    }
    // Add observedKeys to the end of the ordering
    if(!found && observedKeys.find(key) != observedKeys.end()) {
      constrainedGroupBuckets[INT_MAX].push_back(key);
      found = true;
    }
    // Add all other keys to right before the observedKeys
    if(!found) {
      constrainedGroupBuckets[INT_MAX - 1].push_back(key);
    }
  }

  vector<int> cmember(nVars);

  int group = 0;
  for(const auto&[_, bucket] : constrainedGroupBuckets) {
    for(RemappedKey key : bucket) {
      if(key != 0) {
        int index = markedKeysIndex.at(key);
        assert(index < nVars);
        cmember[index] = group;
      }
    }
    group++;
  }

  //double* knobs = nullptr; /* colamd arg 6: parameters (uses defaults if nullptr) */
  double knobs[CCOLAMD_KNOBS];
  ccolamd_set_defaults(knobs);
  knobs[CCOLAMD_DENSE_ROW] = -1;
  knobs[CCOLAMD_DENSE_COL] = -1;

  int stats[CCOLAMD_STATS]; /* colamd arg 7: colamd output statistics and error codes */

  // call colamd, result will be in p
  /* returns (1) if successful, (0) otherwise*/
  if (nVars > 0) {
    int rv = ccolamd((int) nFactors, (int) nVars, (int) Alen, &A[0], &p[0],
        knobs, stats, &cmember[0]);
    if (rv != 1) {
      throw runtime_error("ccolamd failed with return value " + to_string(rv));
    }
  }

  partialOrdering->resize(nVars);
  for(int i = 0; i < p.size(); i++) {
    if(p[i] == -1) {
      break;
    }
    RemappedKey key = markedKeysVector[p[i]];
    (*partialOrdering)[i] = key;

    assert(key != 0);
    assert(cliques_[key]->marked());
  }

  // Manually add back key 0
  partialOrdering->push_back(0);
  assert(partialOrdering->size() == markedKeys.size());

  // cout << "Ordering: ";
  // for(auto k : *partialOrdering) {
  //   cout << (int) unmapKey(k) << " ";
  // }
  // cout << endl;
}

void CholeskyEliminationTree::updateOrdering(const RemappedKeySet& markedKeys, 
                                             const RemappedKeySet& observedKeys,
                                             const ISAM2UpdateParams& params) {
  // cout << "[CholeskyEliminationTree] updateOrdering()" << endl;
  postOrder_ = true;
  orderingVersion_++;
  vector<RemappedKey> partialOrdering;
  getPartialReordering(markedKeys, observedKeys, params, &partialOrdering);
  assert(partialOrdering.size() == markedKeys.size());

  // Adjust ordering_. Shift not reordered keys to the front
  auto minIt = min_element(partialOrdering.begin(), partialOrdering.end(), orderingLess_);
  RemappedKey minReorderedKey = *minIt;
  size_t minReorderedIndex = keyToOrdering_[minReorderedKey];
  size_t fixedIndex = minReorderedIndex;
  vector<RemappedKey> newOrderingToKey(orderingToKey_);
  for(size_t i = minReorderedIndex; i < orderingToKey_.size(); i++) {
    RemappedKey key = orderingToKey_[i];
    if(!cliques_[key]->marked()) {
      newOrderingToKey[fixedIndex] = key;
      fixedIndex++;
    }
  }
  
  // Add newly reordered keys to the back of the ordering
  assert(fixedIndex + partialOrdering.size() == orderingToKey_.size());
  for(size_t i = 0; i < partialOrdering.size(); i++) {
    RemappedKey key = partialOrdering[i];
    newOrderingToKey[fixedIndex + i] = key;
  }
  orderingToKey_ = newOrderingToKey;

  // Now populate keyToOrdering_
  for(size_t i = minReorderedIndex; i < orderingToKey_.size(); i++) {
    RemappedKey key = orderingToKey_[i];
    keyToOrdering_[key] = i;
  }

  unordered_set<sharedFactorWrapper> affectedFactors;

  // Reset all affected nodes
  for(RemappedKey key : partialOrdering) {
    sharedNode node = nodes_[key];
    sharedClique clique = cliques_[key];
    assert(clique->orderingVersion != orderingVersion_);
    clique->orderingVersion = orderingVersion_;

    // If the node has any children cliques, the children need to be reattached
    // Since they might have new parents
    // Need to copy as changing node->clique->children invalidate iterators
    vector<sharedClique> cliqueChildren;
    for(sharedClique child : clique->children) {
      cliqueChildren.push_back(child);
    }
    for(sharedClique child : cliqueChildren) {
      child->reorderAndFindParent();
    }

    for(auto it = node->factors.begin(); it != node->factors.end();) {
      sharedFactorWrapper factor = *it;
      assert(factor != nullptr);
      if(factor->status() == REMOVING || factor->status() == REMOVED) {
        // TODO: Use a better way of removing. Currently need to still be in order
        // because updateOrdering requires factors to be in order
        // Update lambda structure
        node->removeFactor(factor);
        it = node->factors.erase(it);
      }
      else {
        affectedFactors.insert(factor);
        it++;
      }
    }
  }

  // TODO: 
  // Remove all REMOVING factors since all of their keys would have been marked;
  // therefore reodered
  // Remove REMOVED factors from marked keys. There may be REMOVED factors that contain
  // both marked keys and unmarked keys. 
  // No need to reorder factor keys, but need to update the lowest key
  // of a factor
  for(auto factor : affectedFactors) {
    factor->updateLowestHighestKeys();
  }

}

void CholeskyEliminationTree::postReordering(const vector<RemappedKey>& partialOrdering) {
  vector<RemappedKey> newOrderingToKey(orderingToKey_);

  RemappedKey minReorderedKey = partialOrdering.front();
  size_t minReorderedIndex = keyToOrdering_[minReorderedKey];

  assert(keyToOrdering_[0] == keyToOrdering_.size() - 1);
  assert(orderingToKey_.back() == 0);
  assert(partialOrdering.back() == 0);

  size_t curIndex = minReorderedIndex;
  for(size_t i = minReorderedIndex; i < orderingToKey_.size(); i++) {
    RemappedKey key = orderingToKey_[i];
    sharedClique clique = cliques_[key];

    // If come accross the first key of a clique, reorder all keys in the clique
    // to be right behind the first key
    if(key == clique->frontKey()) {
      for(size_t i = 0; i < clique->cliqueSize(); i++) {
        RemappedKey cliqueKey = get<BLOCK_INDEX_KEY>(clique->blockIndices[i]);
        assert(curIndex < newOrderingToKey.size());
        newOrderingToKey[curIndex] = cliqueKey;
        curIndex++;
      }
    }
  }
  assert(curIndex == keyToOrdering_.size());
  orderingToKey_ = newOrderingToKey;

  for(size_t i = minReorderedIndex; i < orderingToKey_.size(); i++) {
    RemappedKey key = orderingToKey_[i];
    keyToOrdering_[key] = i;
  }

}

void CholeskyEliminationTree::allocateStack() {
  if(postOrder_) {
    allocateStackPostOrdering();
  }
  else {
    allocateStackRegular();
  }

  workspace_.allocate(root_->accumSize);
}

void CholeskyEliminationTree::allocateStackRegular() {

  // There should be different behaviors for after reordering and before reordering
  vector<pair<sharedClique, bool>> stack(1, {root_, false});
  while(!stack.empty()) {
    auto& curPair = stack.back();
    sharedClique clique = curPair.first;
    bool& expanded = curPair.second;

    if(!expanded) {

      expanded = true;

      if(clique->marked()) {
        // Defaults to reconstructing
        // TODO: We do not have to default to reconstructing
        clique->setStatusReconstruct();

        size_t diagonalHeight = clique->diagonalHeight();
        size_t subdiagonalHeight = clique->subdiagonalHeight();
        size_t totalHeight = diagonalHeight + subdiagonalHeight;
        // clique->selfSize = diagonalHeight * subdiagonalHeight 
        //                     + subdiagonalHeight * subdiagonalHeight;
        clique->selfSize = totalHeight * totalHeight;

        // expand to children only if marked
        for(sharedClique childClique : clique->children) {
          assert(childClique->status() == UNMARKED);
          stack.push_back({childClique, false});
          childClique->accumSize = 0;
        }
      }
      else {
        assert(clique->status() == UNMARKED);

        // If we reordered this iteration, then set everything to RECONSTRUCT 
        for(size_t i = clique->cliqueSize(); i < clique->blockIndices.size() - 1; i++) {
          // For any marked clique, if any of its decendants is unmarked, set to EDIT
          const auto&[key, row, height] = clique->blockIndices[i];
          assert(nodes_[key]->clique()->marked());
          nodes_[key]->clique()->setStatusEdit();
        }
      }
    }
    else {
      stack.pop_back();

      if(clique->marked()) {
        // Set node's markedStatus to be the same as the clique
        clique->setNodeStatus();
      }
      else {
        assert(clique->status() == UNMARKED);
      }

      // Only count marked children clique and cliques that have reconstruct 
      // for mem allocation
      size_t maxChildSize = 0;
      for(sharedClique childClique : clique->children) {
        if(childClique->accumSize > maxChildSize) {
          maxChildSize = childClique->accumSize;
        }
      }

      clique->accumSize = clique->selfSize + maxChildSize;
    }
  }
}

void CholeskyEliminationTree::allocateStackPostOrdering() {

  vector<pair<sharedClique, bool>> stack(1, {root_, false});
  while(!stack.empty()) {
    auto& curPair = stack.back();
    sharedClique clique = curPair.first;
    bool& expanded = curPair.second;

    if(!expanded) {

      expanded = true;

      if(clique->marked() || clique->hasMarkedAncestor()) {
        if(clique->marked()) {
          clique->setStatusReconstruct();
        }

        // Allocate workspace for marked cliques and cliques that have RECONSTRUCT ancestors
        size_t diagonalHeight = clique->diagonalHeight();
        size_t subdiagonalHeight = clique->subdiagonalHeight();
        size_t totalHeight = diagonalHeight + subdiagonalHeight;

        clique->selfSize = totalHeight * totalHeight;

        // expand to children 
        for(sharedClique childClique : clique->children) {
          assert(childClique->status() == UNMARKED);
          stack.push_back({childClique, false});
          childClique->accumSize = 0;
        }
      }
    }
    else {
      stack.pop_back();

      if(clique->marked()) {
        // Set node's markedStatus to be the same as the clique
        clique->setNodeStatus();
      }
      else {
        assert(clique->status() == UNMARKED);
      }

      // Only count marked children clique and cliques that have reconstruct 
      // for mem allocation
      size_t maxChildSize = 0;
      for(sharedClique childClique : clique->children) {
        if(childClique->accumSize > maxChildSize) {
          maxChildSize = childClique->accumSize;
        }
      }

      clique->accumSize = clique->selfSize + maxChildSize;
    }
  }
}

void CholeskyEliminationTree::deallocateStack() {
  workspace_.deallocate();
}

void CholeskyEliminationTree::allocateAndGatherClique(sharedClique clique, bool allocate, bool reconstruct) {
  // Note: This function is run with the root clique because we want the root
  // clique to allocate space on the workspace, so that mergeWorkspaceClique can work
  // cout << "[CholeskyEliminationTree] allocateAndGatherClique() " << *clique << " reconstruct = " << reconstruct << endl;

  // At this point clique->blockIndices should already be set up
  BlockIndexVector& blockIndices = clique->blockIndices;

  assert(get<BLOCK_INDEX_KEY>(blockIndices.back()) == 0);

  size_t totalHeight = clique->height();
  size_t diagWidth = clique->width();
  size_t subdiagonalHeight = clique->subdiagonalHeight();

  // clique->workspaceIndex = workspace_.push_matrices({{totalHeight, diagWidth}, 
  //                                                    {subdiagonalHeight, subdiagonalHeight}});
  clique->workspaceIndex = workspace_.push_matrices({{totalHeight, totalHeight}});

  if(!reconstruct) {
    // Gather columns from wherever they is living now
    // Need to do a scatter operation based on the old blockIndices on to the new one
    // The clique who owns the data shall be responsible for the blockIndices
    double* matrixData = workspace_.get_ptrs(clique->workspaceIndex)[0];
    CliqueColumns workspaceColumn(matrixData, &blockIndices);

    for(auto& gatherSource : clique->gatherSources) {
      workspaceColumn.addCliqueColumns(gatherSource);
    }
  }

  for(LocalCliqueColumns& gatherSource : clique->gatherSources) {
    if(!gatherSource.ownsData()) {
      gatherSource.deallocateSelfData();
    }
  }

  if(allocate) {
    // Doing this part in scatterClique()
    // clique->gatherSources.clear();
    // auto matrixSource = make_shared<vector<double>>(totalHeight * diagWidth);
    // clique->gatherSources.push_back(LocalCliqueColumns(matrixSource, &blockIndices));
  }
  else {
    assert(clique->gatherSources.size() == 1);
    assert(clique->ownsColumns());
    assert(clique->gatherSources.front().rows() == totalHeight);
    assert(clique->gatherSources.front().cols() == diagWidth);
    assert(clique->gatherSources.front().startIndex() == 0);
  }

}

void CholeskyEliminationTree::choleskyElimination(const Values& theta) {
  // cout << "[CholeskyEliminationTree] choleskyElimination()" << endl;
  
  vector<pair<sharedClique, bool>> stack(1, {root_, false});
  while(!stack.empty()) {
    auto& curPair = stack.back();
    sharedClique clique = curPair.first;
    sharedClique parent = clique->parent();
    bool& expanded = curPair.second;
    if(!expanded) {
      expanded = true;

      // cout << "Restore pass: " << *clique << " markedStatus = " << clique->status() << endl;

      bool expandChildren = clique->marked();

      if(clique->marked()) {
        // check has edit
        vector<RemappedKey> editCols;
        clique->checkEditOrReconstruct(EDIT, &editCols);
        bool gathered = false;
        // if has edit
        if(!editCols.empty()) {
          // cout << "EditCols: ";
          // for(Key k : editCols) {
          //     cout << k << " ";
          // }
          // cout << endl;

          allocateAndGatherClique(clique, true, false);
          gathered = true;
          editOrReconstructFromClique(clique, editCols, 1);
        }
        if(clique->status() == EDIT) {
          if(!gathered) {
            allocateAndGatherClique(clique, true, false);
            gathered = true;
          }
          restoreClique(clique);
        }
        else {
          if(!gathered) {
            allocateAndGatherClique(clique, true, true);
            gathered = true;
          }
          else {
            resetCliqueColumns(clique);
          }
        }
      }
      else {
        // Reorder fixed nodes
        if(clique->orderingVersion != orderingVersion_) {
          // Reorder each node in the clique
          clique->reorderClique();
        }

        vector<Key> reconstructCols;
        clique->checkEditOrReconstruct(RECONSTRUCT, &reconstructCols);

        if(!reconstructCols.empty()) {
          // cout << "ReconstructCols: ";
          // for(Key k : reconstructCols) {
          //     cout << k << " ";
          // }
          // cout << endl;
          allocateAndGatherClique(clique, false, false);
          editOrReconstructFromClique(clique, reconstructCols, -1);

          expandChildren = true;
        }
      }

      // only add child if it's marked or has reconstruct
      // It is possible to be marked but all ancestors are edits
      if(expandChildren) {
        for(sharedClique child : clique->children) {
          stack.push_back({child, false});
        }
      }
      else {
        stack.pop_back();
      }
    }
    else {
      // Eliminate pass
      stack.pop_back();

      // cout << "Eliminate pass: " << *clique << " markedStatus = " << clique->status() << endl;

      if(!clique->marked()) {
        assert(postOrder_ && clique->hasMarkedAncestor());
      }

      // All factors of this clique should be relinearized
      // Do AtA for each node
      // cout << "before construct lambda clique" << endl;
      constructLambdaClique(clique, theta);
      
      if(clique->marked()) {
        // For UNAMRKED nodes, Add AtA blocks for reconstruct columns 
        // that may be connected to unmarked nodes
        // As of Jul. 2023, we do not need to do this for UNMARKED nodes because 
        // we do not depend on UNMARKED nodes

        // Eliminiate clique
        // cout << "before eliminate clique" << endl;
        eliminateClique(clique);
      }

      // Merge with parent
      // cout << "before merge clique" << endl;
      mergeWorkspaceClique(clique);

      // Scatter workspace back into columns
      // UNMARKED nodes still need to scatter to clear workspace
      // cout << "before scatter clique" << endl;
      scatterClique(clique);

      // cout << "before reset clique" << endl;
      clique->resetAfterCholesky();
      // cout << "after reset clique" << endl;

      if(clique->marked()) {
        clique->setBacksolve(true);
      }
      else {
        clique->setBacksolve(false);
      }
    }
  }

  assert(workspace_.empty());

  // checkInvariant_afterCholesky();
}

void CholeskyEliminationTree::editOrReconstructFromClique(
  sharedClique clique,
  const vector<Key>& destCols,
  double sign) {
  if(clique->isLastRow()) { return; }

  // cout << "editOrReconstructFromClique clique: " << *clique << endl;
  // cout << "cols = ";
  // for(Key key : destCols) {
  //   cout << key << " ";
  // }
  // cout << endl;

  auto m = workspace_.get_matrices(clique->workspaceIndex)[0];

  size_t cliqueSize = clique->cliqueSize();
  const auto& blockIndices = clique->blockIndices;
  size_t totalHeight = clique->height();
  size_t diagWidth = clique->width();
  size_t subdiagHeight = clique->subdiagonalHeight();

  assert(m.rows() == totalHeight && m.cols() == totalHeight);

  auto B = Eigen::Block<Eigen::Map<ColMajorMatrix>>
            (m, diagWidth, 0, subdiagHeight, diagWidth);
  auto C = Eigen::Block<Eigen::Map<ColMajorMatrix>>
            (m, diagWidth, diagWidth, subdiagHeight, subdiagHeight);

  assert(get<BLOCK_INDEX_KEY>(blockIndices.back()) == 0);

  bool processGrouped = true;
  size_t destSize = destCols.size();
  size_t totalSize = blockIndices.size() - 1;
  // if(destSize * 2 < totalSize) {
  //     processGrouped = false;
  // }

  if(processGrouped) {

    C.selfadjointView<Eigen::Lower>().rankUpdate(B, sign);

    auto destIt = destCols.begin();
    auto destEnd = destCols.end();

    // Set unused columns to be 0 since we did everything at once
    for(size_t i = cliqueSize; i < blockIndices.size(); i++) {
      const auto&[key, col, width] = blockIndices[i];

      if(destIt == destEnd || *destIt != key) {
        block(m, col, col, totalHeight - col, width).setZero();
      }
      else if (destIt != destEnd){
        destIt++;
      }
    }

    assert(destIt == destEnd);
  }
  else {
    auto destIt = destCols.begin();
    auto destEnd = destCols.end();

    for(size_t i = cliqueSize; i < blockIndices.size(); i++) {
      RemappedKey key = get<BLOCK_INDEX_KEY>(blockIndices[i]);   

      if(destIt != destEnd && *destIt == key) {
        destIt++;
        size_t col = get<BLOCK_INDEX_ROW>(blockIndices[i]);
        size_t width = get<BLOCK_INDEX_HEIGHT>(blockIndices[i]);
        size_t cHeight = totalHeight - col;
        block(m, col, col, cHeight, width) += 
          sign * block(m, col, 0, cHeight, diagWidth)
          * block(m, col, 0, width, diagWidth).transpose();
      }
    }
  }
}

void CholeskyEliminationTree::restoreClique(sharedClique clique) {
  if(clique->isLastRow()) { return; }

  auto m = workspace_.get_matrices(clique->workspaceIndex)[0];

  size_t totalHeight = clique->height();
  size_t diagWidth = clique->width();
  auto D = block(m, 0, 0, diagWidth, diagWidth);
  auto DB = block(m, 0, 0, totalHeight, diagWidth);

  D.triangularView<Eigen::StrictlyUpper>().setZero();

  DB *= D.transpose();
}

void CholeskyEliminationTree::resetCliqueColumns(sharedClique clique) {
  if(clique->isLastRow()) { return; }

  Eigen::Map<ColMajorMatrix> m = workspace_.get_matrices(clique->workspaceIndex)[0];

  size_t totalHeight = clique->height();
  size_t diagWidth = clique->width();

  auto DB = Eigen::Block<Eigen::Map<ColMajorMatrix>>(m, 0, 0, totalHeight, diagWidth);
  DB.setZero();
}

struct FactorInfo {
  // <row, ordering, status>
  bool valid = false;
  vector<tuple<size_t, size_t, MarkedStatus>> v;
  void construct(
      const vector<RemappedKey>& factorKeys,
      const unordered_map<RemappedKey, tuple<size_t, size_t, MarkedStatus>>& infoMap) {
    valid = true;
    v.reserve(factorKeys.size());
    for(const auto& key : factorKeys) {
      v.push_back(infoMap.at(key));
    }
  }
  const tuple<size_t, size_t, MarkedStatus>& operator()(const size_t& i) const {
    return v[i];
  }
};

struct CheckEdits {
  bool operator()(const MarkedStatus& status) const { 
    return status == EDIT;
  }
};
struct CheckNotUnmarked {
  bool operator()(const MarkedStatus& status) const { 
    return !(status == UNMARKED);
  }
};
struct CheckReconstructNew {
  bool operator()(const MarkedStatus& status) const { 
    return (status == RECONSTRUCT) || (status == NEW);
  }
};

void CholeskyEliminationTree::constructLambdaClique(sharedClique clique, const Values& theta) {
  if(clique->isLastRow()) { return; }

  // cout << "constructLambda for clique " << *clique << " status = " << clique->status() << endl;

  auto m = workspace_.get_matrices(clique->workspaceIndex)[0];

  const auto& blockIndices = clique->blockIndices;

  unordered_map<RemappedKey, tuple<size_t, size_t, MarkedStatus>> infoMap;
  for(const auto&[key, row, height] : blockIndices) {
    assert(nodes_[key]->status() != MARGINALIZED);
    infoMap.insert({key, {row, keyToOrdering_[key], nodes_[key]->status()}});
  }

  for(sharedNode node : clique->nodes) {

    // There are 4 types of MarkedStatus and 6 types of FactorStatus
    // enum MarkedStatus {UNMARKED, EDIT, RECONSTRUCT, NEW};
    // enum FactorStatus {UNLINEARIZED, LINEARIZED, RELINEARIZE, LINEAR, REMOVING, REMOVED};
    // We only process a factor if the lowest key of that factor is in our clique
    // so that memory for that factor is guaranteed to be allocated
    // The MarkedStatus of each node in a clique may be different, generally ony happening 
    // if a clique is EDIT and we add some NEW nodes at the end
    // The MarkedStatus of each node in a factor may also be different, because
    // the nodes may be in different cliques
    // 1. Subtract all MarkedStatus-FactorStatus that are EDIT-{RELINEARIZE, REMOVING}
    // 2. Linearize/Relienarize all factors that are {UNLINEARIZED, RELINEARIZE}
    // 3. Add all MarkedStatus-FactorStatus that are {EDIT, RECONSTRUCT, NEW}-{UNLINEARIZED, RELINEARIZE} + {RECONSTRUCT, NEW}-{LINEARIZED, LINEAR}

    for(sharedFactorWrapper factorWrapper : node->factors) {
      // Node must be lowest node in factor
      // We need this to ensure there is an allocated entry in the workspace
      // for this AtA bock
      if(node->key != factorWrapper->lowestKey()) {
        continue;
      }
      else if(factorWrapper->status() == REMOVED) {
        continue;
      }

      // Status of the factor before constructing lambda
      // This allows us to change factorWrapper->status() without worrying about the logic
      FactorStatus factorStatus = factorWrapper->status();

      FactorInfo info;

      // Handle edits first
      // Subtract Hessian of the cachedLinearFactor from the workspace
      if(factorStatus == RELINEARIZE || factorStatus == REMOVING) {
        sharedLinearFactor linearFactor = factorWrapper->getLinearFactor();
        assert(linearFactor != nullptr);

        info.construct(factorWrapper->remappedKeys(), infoMap);

        factorWrapper->updateHessian(m, -1, info, CheckEdits());
      }

      if(factorStatus == REMOVING) { 
        factorWrapper->setStatusRemoved(); 
        continue;
      }

      // Then relinearize factor to linear factor
      bool unlinearized = factorWrapper->linearizeIfNeeded(theta);

      // Then add hessian of factor to workspace
      if(unlinearized) {
        // // DEBUG
        // cout << "ADD factor " << factorWrapper << " keys: ";
        // factorWrapper->printKeys(cout);
        // cout << endl;
        // // END DEBUG
        if(!info.valid) {
          info.construct(factorWrapper->remappedKeys(), infoMap);
        }
        factorWrapper->updateHessian(m, 1, info, CheckNotUnmarked());
      }
      else if(factorStatus == LINEARIZED || factorStatus == LINEAR) {
        // // DEBUG
        // cout << "ADD factor " << factorWrapper << " keys: ";
        // factorWrapper->printKeys(cout);
        // cout << endl;
        // // END DEBUG
        if(!info.valid) {
          info.construct(factorWrapper->remappedKeys(), infoMap);
        }
        factorWrapper->updateHessian(m, 1, info, CheckReconstructNew());
      }
    }
  }
}

void CholeskyEliminationTree::eliminateClique(sharedClique clique) {
  if(clique->isLastRow()) { return; }

  auto m = workspace_.get_matrices(clique->workspaceIndex)[0];

  size_t totalHeight = clique->height();
  size_t bWidth = clique->width();
  size_t bHeight = totalHeight - bWidth;

  assert(totalHeight == m.rows());

  // // DEBUG
  // cout << "Clique: " << *clique << endl << endl;
  // m.triangularView<Eigen::StrictlyUpper>().setZero();
  // cout << "Before eliminate. m = \n" << block(m, 0, 0, totalHeight, bWidth) << endl << endl;

  Eigen::Block<Eigen::Map<ColMajorMatrix>> D = block(m, 0, 0, bWidth, bWidth);
  Eigen::LLT<Eigen::Ref<Eigen::Block<Eigen::Map<ColMajorMatrix>>>> llt(D);
  if(llt.info() == Eigen::NumericalIssue) {
    cout << "Diagonal block not positive definite!" << endl;
    clique->printClique(cout);
    exit(1);
  }
  auto L = D.triangularView<Eigen::Lower>();
  auto B = block(m, bWidth, 0, bHeight, bWidth);
  auto C = block(m, bWidth, bWidth, bHeight, bHeight);

  L.solveInPlace(B.transpose());

  if(bHeight != -1) {
    C.selfadjointView<Eigen::Lower>().rankUpdate(B, -1);
  }

  // // DEBUG
  // m.triangularView<Eigen::StrictlyUpper>().setZero();
  // cout << "After eliminate. m = \n" << block(m, 0, 0, totalHeight, bWidth) << endl << endl;
}

void CholeskyEliminationTree::mergeWorkspaceClique(sharedClique clique) {
  if(clique->isLastRow()) { 
    assert(clique == root_);
    assert(!clique->parent());
    return; 
  }

  sharedClique parent = clique->parent();

  // Interpret the workspace matrices as a CliqueColumns, then merge them
  auto childMatrixData = workspace_.get_ptrs(clique->workspaceIndex)[0];
  CliqueColumns childCliqueColumns(childMatrixData, 
      &clique->blockIndices, 
      clique->cliqueSize());

  auto parentMatrixData = workspace_.get_ptrs(parent->workspaceIndex)[0];
  CliqueColumns parentCliqueColumns(parentMatrixData, &parent->blockIndices);

  parentCliqueColumns.addCliqueColumns(childCliqueColumns);
}

void CholeskyEliminationTree::scatterClique(sharedClique clique) {

  if(clique->marked()) {
    auto m_ptr = workspace_.get_ptrs(clique->workspaceIndex)[0];
    size_t r = clique->height(), c = clique->width();
    std::shared_ptr<vector<double>> matrixSource = std::make_shared<vector<double>>(r * c);
    std::shared_ptr<BlockIndexVector> blockIndicesSource 
      = std::make_shared<BlockIndexVector>(clique->blockIndices);

    // // DEBUG
    // auto m = workspace_.get_matrices(clique->workspaceIndex)[0];
    // cout << "m rows = " << m.rows() << " m cols = " << m.cols() << " r = " << r << " c = " << c << endl;
    // assert(r == m.rows() && c <= m.cols());

    memcpy(matrixSource->data(), m_ptr, sizeof(double) * r * c);

    clique->gatherSources.clear();
    clique->gatherSources.push_back(LocalCliqueColumns(matrixSource, 
                                                       blockIndicesSource, 
                                                       0, clique->cliqueSize()));
  }

  workspace_.pop_matrices();
}

void CholeskyEliminationTree::backsolve(VectorValues* delta_ptr, double tol) {
  // cout << "[CholeskyEliminationTree] backwardSolve()" << endl;
  // Do a pre-order traversal from top ot bottom
  // For each node, first process the belowDiagonalBlocks, then do solve on the transpose of the diagonal
  vector<pair<sharedClique, bool>> stack(1, {root_, false});
  // Already solved delta. Used to avoid VectorValues.at(key)
  vector<Vector> cachedDelta(nodes_.size());
  while(!stack.empty()) {
    auto& curPair = stack.back();
    sharedClique clique = curPair.first;
    bool& expanded = curPair.second;

    if(!expanded) {
      expanded = true;
      if(clique->orderingVersion != orderingVersion_) {
        clique->reorderClique();
      }

      bool backsolve = clique->marked() || clique->needsBacksolve();

      if(backsolve) {

        backsolveClique(clique, delta_ptr, tol);

        for(sharedClique child : clique->children) {
          stack.push_back({child, false});
        }
      }
    }
    else {
      stack.pop_back();

      clique->resetAfterBacksolve();
    }
  }

  deallocateStack();

#ifdef DEBUG
  checkInvariant_afterBackSolve();
#endif
}

void CholeskyEliminationTree::backsolveClique(
    sharedClique clique, 
    VectorValues* delta_ptr, 
    double tol) {
  if(clique->isLastRow()) { return; }

  // cout << "CholeskyEliminationTree::backsolveClique(): " << *clique << endl;

  // Columns are already gathered
  assert(clique->gatherSources.size() == 1);

  const auto& blockIndices = clique->blockIndices;
  auto m = clique->gatherSources.front().matrix();
  size_t cliqueSize = clique->cliqueSize();
  size_t totalHeight = clique->height();
  size_t diagWidth = clique->width();
  size_t subdiagHeight = clique->subdiagonalHeight();

  // Copy over L^-1 Atb row into delta
  Vector delta = block(m, totalHeight - 1, 0, 1, diagWidth).transpose();

  if(subdiagHeight > 1) {
    Vector gatherX(subdiagHeight - 1);

    // Gather deltas this clique depends on, don't need last row
    for(int i = cliqueSize; i < blockIndices.size() - 1; i++) {
      auto[key, row, height] = blockIndices[i];
      Key unmappedKey = unmapKey(key);

      row -= diagWidth;

      gatherX.block(row, 0, height, 1) = delta_ptr->at(unmappedKey);
    }

    auto B = block(m, diagWidth, 0, subdiagHeight - 1, diagWidth); // sub-diagonal blocks
    delta -= B.transpose() * gatherX;
  }
  
  // Solve diagonal
  auto D = block(m, 0, 0, diagWidth, diagWidth);
  auto LT = D.transpose().triangularView<Eigen::Upper>();
  LT.solveInPlace(delta);

  // Copy delta back into delta_ptr
  for(int i = 0; i < cliqueSize; i++) {
    auto[key, row, height] = blockIndices[i];
    Key unmappedKey = unmapKey(key);

    Vector delta_diff = delta.block(row, 0, height, 1) - delta_ptr->at(unmappedKey);
    delta_ptr->at(unmappedKey) = delta.block(row, 0, height, 1);

    if(valuesChanged(delta_diff, tol)) {
      nodes_.at(key)->backsolve = true;
    }
  }
}

bool CholeskyEliminationTree::valuesChanged(const Vector& diff, double tol) {
  return diff.lpNorm<Eigen::Infinity>() >= tol;
}

void CholeskyEliminationTree::reset() {
  postOrder_ = false;
}

// FIXME: do this for a set of keys, not just a single key
void CholeskyEliminationTree::getAffectedKeys(Key unmappedKey, set<Key>& additionalKeys) const {
  // cout << "in getAffectedKeys, cliques_.size() = " << cliques_.size() << endl;

  RemappedKey startKey = getRemapKey(unmappedKey);
  sharedClique curClique = cliques_[startKey];

  // Everything above index should be affected
  size_t index = -1;
  for(index = 0; index < curClique->cliqueSize(); index++) {
    RemappedKey key = get<BLOCK_INDEX_KEY>(curClique->blockIndices[index]);
    if(key == startKey) { break; }
  }
  assert(index != -1);

  do {
    for(; index < curClique->cliqueSize(); index++) {
      RemappedKey key = get<BLOCK_INDEX_KEY>(curClique->blockIndices[index]);
      Key unmappedKey = unmapKey(key);
      additionalKeys.insert(unmappedKey);
    }

    // Start from index 0 in parent
    curClique = curClique->parent();
    index = 0;
  } while(curClique != root_);
}

void CholeskyEliminationTree::getAffectedDescendantKeys(
    Key unmappedKey, std::set<Key>& additionalKeys) const {
  
  RemappedKey key = getRemapKey(unmappedKey);

  vector<sharedClique> stack({cliques_[key]});

  while(!stack.empty()) {
    sharedClique curClique = stack.back();
    stack.pop_back();

    for(sharedNode node : curClique->nodes) {
      additionalKeys.insert(unmapKey(node->key));
    }

    // TODO: optimize this. If additionalKeys already has key, then we don't need to check
    // descendants

    for(sharedClique child : curClique->children) {
      stack.push_back(child);
    }
  }
}

void CholeskyEliminationTree::marginalizeLeaves(
    const FastList<Key>& leafKeys,
    boost::optional<FactorIndices&> marginalFactorsIndices,
    boost::optional<FactorIndices&> deletedFactorsIndices) {

  // First remap keys and mark all nodes we want to marginalize
  vector<sharedClique> marginalizedCliques;
  vector<RemappedKey> remappedLeafKeys;
  for(Key unmappedKey : leafKeys) {
    RemappedKey key = getRemapKey(unmappedKey);
    remappedLeafKeys.push_back(key);
    assert(nodes_[key]->status() == UNMARKED);
    nodes_[key]->setStatusMarginalized();
    sharedClique clique = cliques_[key];
    if(key == clique->frontKey()) {
      // No need to worry about duplicates if we only mark clique 
      // for the first key of the clique. If a later key in the clique is marginalized
      // the first key must also be marginalized
      clique->setStatusMarginalized();
      marginalizedCliques.push_back(clique);
    }
  }

  // Remove cliques that have no connections to any unmarginalized nodes
  for(int i = marginalizedCliques.size() - 1; i >= 0; i--) {
    sharedClique clique = marginalizedCliques[i];
    RemappedKey lastKey = get<BLOCK_INDEX_KEY>(*(clique->blockIndices.rbegin() + 1));
    if(nodes_[lastKey]->status() == MARGINALIZED) {
      // clique can be safely removed
      marginalizedCliques[i] = marginalizedCliques.back();
      marginalizedCliques.pop_back();
    }
  }

  // Mark all the nodes that we want to marginalize
  FactorIndices deletedFactors;
  FactorIndices marginalFactors;
  for(sharedClique clique : marginalizedCliques) {
    marginalizeClique(clique, &deletedFactors, &marginalFactors);
  }


#ifdef DEBUG
  checkInvariant_afterMarginalize();
#endif
}

void CholeskyEliminationTree::marginalizeClique(
      sharedClique clique,
      FactorIndices* deletedFactors,
      FactorIndices* marginalFactors) {
  // cout << "[CholeskyEliminationTree] marginalizedClique() clique: " << *clique << endl;

  // Set all factors in each node to be LINEAR
  // and remove all keys that are marginalized out
  // If all the keys in a factor is marginalized out, remove that factor
  // Remove all factors in each node
  // There is a catch here. If a factor is set to REMOVING and its lower keys 
  // are marginalized out, we need to update the factor so that its lowest key
  // is the lowest remaining key. This is because if that key is set to EDIT
  // later on, the factor needs to be subtracted. We can simple remove the block indices
  // corresponding to the marginalized keys
  for(sharedNode node : clique->nodes) {
    if(node->status() != MARGINALIZED) { continue; }
    for(sharedFactorWrapper factor : node->factors) {
      if(factor->lowestKey() == node->key) {
        bool deleteFactor = factor->marginalizeKeys();
        if(deleteFactor) {
          deletedFactors->push_back(factor->factorIndex());
        }
        else {
          assert(factor->status() == LINEAR);
        }
      }
    }
  }

  assert(clique->gatherSources.size() == 1 && clique->ownsColumns());
  LocalCliqueColumns cliqueColumns = clique->gatherSources.front();

  const BlockIndexVector& blockIndices = clique->blockIndices;
  size_t lowestIndex = 0, lowestRow = 0;
  for(size_t i = 0; i < blockIndices.size() - 1; i++) {
    const auto&[key, row, height] = blockIndices[i];
    if(nodes_[key]->status() != MARGINALIZED) {
      lowestIndex = i;
      lowestRow = row;
      break;
    }
  }
  assert(lowestIndex != 0);

  
  bool splitClique = lowestIndex < clique->cliqueSize();
  size_t splitIndex = lowestIndex;
  size_t totalWidth = clique->width();
  size_t lowestCol = splitClique? lowestRow : totalWidth;
  size_t newHeight = clique->height() - lowestRow;

  auto origCol = cliqueColumns.matrix();
  auto B = block(origCol, lowestRow, 0, newHeight, lowestCol);

  // Construct SymmetricBlockMatrix for HessianFactor
  vector<size_t> dimensions;
  vector<Key> unmappedKeys;
  dimensions.reserve(blockIndices.size());
  for(size_t i = lowestIndex; i < blockIndices.size() - 1; i++) {
    // No need to add -1 key
    const auto&[key, row, height] = blockIndices[i];
    dimensions.push_back(height);
    unmappedKeys.push_back(unmapKey(key));
  }
  SymmetricBlockMatrix factorMatrix(dimensions, true);
  factorMatrix.matrix().setZero();
  assert(B.rows() == factorMatrix.selfadjointView().rows());
  factorMatrix.selfadjointView().rankUpdate(B, -1);

  // Make HessianFactor
  HessianFactor::shared_ptr marginalFactor = boost::make_shared<HessianFactor>(unmappedKeys, factorMatrix);

  size_t factorIndex = factors_.size();
  sharedFactorWrapper factorWrapper = std::make_shared<FactorWrapper>(
                                        factorIndex, nullptr, marginalFactor, this);

  factors_.push_back(factorWrapper);

  for(RemappedKey key : factorWrapper->remappedKeys()) {
    nodes_[key]->addFactor(factorWrapper);
  }

  if(splitClique) {
    // If a clique is split in the middle, We need to create a new Clique with 
    // correct CliqueColumns, and that new clique needs to own its CliqueColumns
    sharedClique newClique = clique->splitClique(splitIndex);
  }

  // Detach clique from parent; mark each node in the clique as MARGINALIZED (fixed). If fixed nodes are relinearized, ignore. Also ignore new factors connected to fixed nodes
  clique->detachParent();

  for(sharedNode node : clique->nodes) {
    if(node->status() != MARGINALIZED) {
      cout << "Clique " << *clique << " nodes are not all MARGINALIZED" << endl;
      for(sharedNode node : clique->nodes) {
        cout << node->key << " " << node->status() << endl;
      }
    }
    assert(node->status() == MARGINALIZED);
  }
}

void CholeskyEliminationTree::printOrderingUnmapped(std::ostream& os) const {
  os << "Unmapped Ordering: ";
  for(RemappedKey key : orderingToKey_) {
    os << (Key) unmapKey(key) << " ";
  }
  os << endl;
}

void CholeskyEliminationTree::printOrderingRemapped(std::ostream& os) const {
  os << "Remapped Ordering: ";
  for(RemappedKey key : orderingToKey_) {
    os << key << " ";
  }
  os << endl;
}

NonlinearFactor::shared_ptr CholeskyEliminationTree::nonlinearFactorAt(FactorIndex factorIndex) { 
  size_t realFactorIndex = factorIndexTransformMap_[factorIndex];
  return factors_[realFactorIndex]->nonlinearFactor(); 
}

RemappedKey CholeskyEliminationTree::addRemapKey(const Key unmappedKey) {

  // Remap regular keys to start from 1, 0 is for last row (unmappedKey -1)
  RemappedKey remappedKey = keyTransformMap_.size();

  bool inserted = keyTransformMap_.insert({unmappedKey, remappedKey}).second;

  if(!inserted) {
    return keyTransformMap_.at(unmappedKey);
  }

  assert(remappedKey == keyToOrdering_.size());

  if(unmappedKey != -1) {
    // Insert the new regular key in front of the last row
    keyToOrdering_[0] = keyToOrdering_.size();
    keyToOrdering_.push_back(keyToOrdering_.size() - 1);
    orderingToKey_.back() = remappedKey;
    orderingToKey_.push_back(0);
  }
  else {
    keyToOrdering_.push_back(0);
    orderingToKey_.push_back(0);
  }

  unmappedKeys_.push_back(unmappedKey);

  return remappedKey;
}

RemappedKey CholeskyEliminationTree::getRemapKey(const Key unmappedKey) const {
  assert(keyTransformMap_.find(unmappedKey) != keyTransformMap_.end());
  return keyTransformMap_.at(unmappedKey);
}

Key CholeskyEliminationTree::unmapKey(const RemappedKey remappedKey) const {
  return unmappedKeys_.at(remappedKey);
}

void CholeskyEliminationTree::addNewNode(const RemappedKey key, const size_t width) {
  sharedNode newNode = make_shared<Node>(this, key, width);
  nodes_.push_back(newNode);
  cliques_.push_back(nullptr);

  // We cannot make sharedNode part of the constructor
  // Because at construct time, Clique is not pointed to by shared_ptr
  // Maybe there is a cleaner way?
  sharedClique newClique = make_shared<Clique>(this);
  newClique->addNode(newNode);

  assert(newNode->clique() == newClique);

  totalDeltaDim_ += width;
}

size_t CholeskyEliminationTree::colWidth(const RemappedKey key) const {
  return nodes_.at(key)->width;
}

void CholeskyEliminationTree::checkInvariant_afterSymbolic() const {
  for(sharedClique clique : cliques_) {
    clique->checkInvariant();
  }

  assert(orderingToKey_.back() == 0);
  assert(keyToOrdering_[0] == keyToOrdering_.size() - 1);

  for(int i = 0; i < orderingToKey_.size(); i++) {
    auto key = orderingToKey_[i];
    assert(keyToOrdering_[key] == i);
  }

}

void CholeskyEliminationTree::checkInvariant_afterBackSolve() const {
  for(const sharedClique& clique : cliques_) {
    if(clique->isLastRow()) {
      assert(clique->parent() == nullptr);
      assert(clique == root_);
      assert(clique.use_count() == 2);  // one from cliques_ and one from root_
    }
    else {
      // nodes.size() from cliques_ and one of clique->children
      assert(clique.use_count() == clique->cliqueSize() + 1);  
    }
    assert(clique->ownsColumns());
    clique->checkInvariant();
  }
  for(const sharedNode& node : nodes_) {
    // One in nodes_ and one in clique->nodes
    assert(node.use_count() == 2);
  }
}

void CholeskyEliminationTree::checkInvariant_afterMarginalize() const {
  cout << "check invariant after marginalized" << endl;
  {
    for(sharedNode node : nodes_) {
      sharedClique clique = node->clique();
      bool flag = false;
      for(sharedNode cliqueNode : clique->nodes) {
        if(cliqueNode == node) {
          flag = true;
          break;
        }
      }
      assert(flag);
    }
  }
  {
    vector<sharedClique> stack(1, root_);

    while(!stack.empty()) {
      sharedClique curClique = stack.back();
      stack.pop_back();

      cout << "checking clique " << *curClique << endl;

      for(sharedClique child : curClique->children) {
        assert(child->status() == UNMARKED);
        stack.push_back(child);
      }
    }
  }
  {
    for(sharedNode node : nodes_) {
      RemappedKey key = node->key;
      assert(node->status() == UNMARKED || node->status() == MARGINALIZED);
      if(node->status() == UNMARKED) {
        sharedClique clique = node->clique();
        while(clique) {
          if(clique->status() != UNMARKED) {
            cout << "node " << node->key << " is UNMARKED, but clique " << *clique << " is status " << clique->status() << endl;
          } 
          assert(clique->status() == UNMARKED);
          sharedClique parent = clique->parent();
          if(parent == nullptr) {
            assert(clique == root_);
          }
          clique = parent;
        }
      }
      else {
        sharedClique clique = node->clique();
        while(clique) {
          assert(clique->status() == MARGINALIZED);
          sharedClique parent = clique->parent();
          if(parent == nullptr) {
            assert(clique != root_);
          }
          clique = parent;
        }
      }
    }
  }
}

} // namespace gtsam
