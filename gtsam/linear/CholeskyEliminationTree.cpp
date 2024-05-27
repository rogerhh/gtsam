
/**
* @file    CholeskyEliminationTree.h
* @brief   Elimination tree structure to perform symbolic factorization and Cholesky factorization
* @author  Roger Hsiao (rogerhh)
* @date    Feb. 8, 2023
*/

#include "gtsam/inference/graph-inl.h"
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
#include <locale>
#include <stdexcept>
#include <cstdlib>
#include <queue>

using namespace std;

namespace gtsam {

// Convenience function
inline Eigen::Block<Eigen::Map<ColMajorMatrix<GEMMINI_TYPE>>> block(
    Eigen::Map<ColMajorMatrix<GEMMINI_TYPE>> m, 
    size_t r, size_t c, size_t h, size_t w) {
    return Eigen::Block<Eigen::Map<ColMajorMatrix<GEMMINI_TYPE>>>(m, r, c, h, w);
}

CholeskyEliminationTree::CholeskyEliminationTree() : orderingLess_(this) {
  // Add new node for the last row. Technically should be connected to everything
  RemappedKey key = addRemapKey(-1);
  assert(key == 0);
  addNewNode(key, 1); 
}

void CholeskyEliminationTree::addVariables(const Values& newTheta) {
  // cout << "[CholeskyEliminationTree] addVariables() " << newTheta.size() << endl;
  datasetgen_newKeys.clear();
  for(const auto& keyValPair : newTheta) {
    const Key& unmappedKey = keyValPair.key;
    const Value& val = keyValPair.value;
    const size_t dim = val.dim();

    // First remap the key to start from 1
    RemappedKey key = addRemapKey(unmappedKey);
    addNewNode(key, dim);

    datasetgen_newKeys.push_back(key);
    
    // Add regularization factor because problem is too ill conditioned
    double lambda = 1e-2; // This needs to be a large number since the Hessian factor takes in the covariance matrix
    Vector mu = Matrix::Zero(dim, 1);
    Matrix sigma = lambda * Matrix::Identity(dim, dim);
    JacobianFactor::shared_ptr regularization_factor 
        = boost::make_shared<JacobianFactor>(unmappedKey, sigma, mu);

    size_t factorIndex = factors_.size();
    sharedFactorWrapper factorWrapper = std::make_shared<FactorWrapper>(
            factorIndex, -1, nullptr, regularization_factor, this);

    factors_.push_back(factorWrapper);

    nodes_[key]->addFactor(factorWrapper);

  }

  // Set this for staleness check in LRU
  for(const auto& keyValPair : newTheta) {
    const Key& unmappedKey = keyValPair.key;
    RemappedKey remappedKey = getRemapKey(unmappedKey);
    nodes_[remappedKey]->lastRelinStep = keyTransformMap_.size();
  }
} 

#define FACTOR7_3_COST 15000
#define FACTOR13_6_COST 20000

int64_t CholeskyEliminationTree::computeCostRelin(
    const RemappedKey remappedKey,
    int num_threads,
    const vector<CostStatus>& curRelinCostStatus) {

  sharedNode node = nodes_[remappedKey];

  if(!node) {
    cout << "In computeCostRelin(), node " << remappedKey << " is null" << endl;
    exit(1);
  }
  assert(node);

  int64_t cost = 0;
  for(sharedFactorWrapper factorWrapper : node->factors) {
    if(factorWrapper->status() == LINEAR) { continue; }
    if(factorWrapper->status() != LINEARIZED) {
      cout << "factor " << *factorWrapper << " status is " << factorWrapper->status() << endl;
      exit(1);
    }
    assert(factorWrapper->status() == LINEARIZED);

    bool addCostFlag = true;
    for(int i = 0; i < factorWrapper->remappedKeys().size(); i++) {
      RemappedKey k = factorWrapper->remappedKeys()[i];
      if(curRelinCostStatus[k] == COST_RELIN) {
        addCostFlag = false;
        break;
      }
    }

    if(addCostFlag) {
      int lastIndex = factorWrapper->blockIndices().size() - 1;
      if(get<BLOCK_INDEX_ROW>(factorWrapper->blockIndices()[lastIndex]) == 12) {
        cost += FACTOR13_6_COST;
      } 
      else {
        cost += FACTOR7_3_COST;
      }
    }
  }
  return cost / num_threads;
}

int64_t CholeskyEliminationTree::computeCost(
    const RemappedKey remappedKey, 
    const int num_threads,
    vector<sharedClique>* updatedCliques) {
  
  // computeCost(key): (assuming reordering)
  // 1. First find the lowest-ordered key connected to the key
  // 2. Go up the tree from the clique
  //     If clique is marked, stop
  //     If clique is fixed, upgrade to marked and compute the cost difference
  //     If clique is unmarked, upgrade to marked and compute cost
  //
  //     Go down the tree from each marked clique,
  //     If descendant is marked, stop (this is an error)
  //     If descendant is fixed, stop
  //     If descendant is unmarked, check if need to be fixed and compute cost if fixed
  //       and propagate to children
  //
  //     while clique is not marked, mark the clique and predict cycles
  //     For each marked clique, find all fixed cliques and predict fixed cost
  //
  //     Parallel case: if clique is parallelizable, reduce cost by 1/overhead_const * p (p is the number of threads) (mplier is set to 0.7 for now)
  //     Root clique is the clique that has root_ as its parent

  sharedClique clique = cliques_[remappedKey];


  int64_t cost = 0;

  while(clique != nullptr && clique != root_) {
    if(clique->costStatus == COST_MARKED || clique->nextCostStatus == COST_MARKED) {
      break;
    }

    clique->computeCostMarked(num_threads);

    if(clique->nextCostStatus == COST_FIXED) { 
      // This check if probably redundant
      int64_t additional_cost = clique->markedCost - clique->fixedCost;
      cost += additional_cost;
    }
    else if(clique->nextCostStatus == COST_UNMARKED) {
      cost += clique->markedCost;
    }

    clique->nextCostStatus = COST_MARKED;
    updatedCliques->push_back(clique);

    // Go down unmarked descendants of the marked node
    vector<sharedClique> stack;
    for(sharedClique childClique : clique->children) {
      if(childClique->nextCostStatus == COST_UNMARKED) {
        stack.push_back(childClique);
      }
    }

    while(!stack.empty()) {
      sharedClique curClique = stack.back();
      stack.pop_back();

      if(curClique->nextCostStatus == COST_UNMARKED) {
        // First check if need to be fixed
        // If second to last key block indices is in this clique
        // we can just check order(secondToLastKey) <= order(clique first key)
        RemappedKey secondToLastKey = get<BLOCK_INDEX_KEY>(curClique->blockIndices[curClique->blockIndices.size() - 2]);
        if(!orderingLess_(secondToLastKey, clique->frontKey())) {
          cost += curClique->computeCostFixed(num_threads);
          curClique->nextCostStatus = COST_FIXED;
          updatedCliques->push_back(curClique);

          stack.insert(stack.end(), curClique->children.begin(), curClique->children.end());
        }
      }
      else if(curClique->nextCostStatus == COST_FIXED) {
        int64_t oldFixedCost = curClique->fixedCost;
        cost += curClique->computeCostFixed(num_threads) - oldFixedCost;
        stack.insert(stack.end(), curClique->children.begin(), curClique->children.end());
      }
      else if(curClique->nextCostStatus == COST_MARKED) {
        cout << "Unexpected cost status for clique: " << *curClique << " status: " << curClique->nextCostStatus << " parent: " << *curClique->parent() << endl;
        exit(1);
      }

    }
    clique = clique->parent();
  }

  return cost;
}

void CholeskyEliminationTree::commitCost(vector<sharedClique>& updatedCliques, 
                                         vector<sharedClique>* allUpdatedCliques) {
  for(sharedClique clique : updatedCliques) {
    clique->costStatus = clique->nextCostStatus;
    allUpdatedCliques->push_back(clique);
  }
}

void CholeskyEliminationTree::uncommitCost(vector<sharedClique>& updatedCliques, 
                                           vector<sharedClique>* allUpdatedCliques) {
  for(sharedClique clique : updatedCliques) {
    clique->nextCostStatus = clique->costStatus;
    allUpdatedCliques->push_back(clique);
  }
}

void CholeskyEliminationTree::resetCost(vector<sharedClique>& updatedCliques) {
  for(sharedClique clique : updatedCliques) {
    clique->resetCost();
  }
}

void CholeskyEliminationTree::checkCliqueParallelizable() {
  sharedClique curClique = root_;
  while(curClique) {
    curClique->parallelizable = false;
    if(curClique->children.size() == 1) {
      curClique = *(curClique->children.begin());
    }
    else {
      break;
    }
  }
}

void CholeskyEliminationTree::resetCliqueParallelizable() {
  sharedClique curClique = root_;
  while(curClique) {
    curClique->parallelizable = true;
    if(curClique->children.size() == 1) {
      curClique = *(curClique->children.begin());
    }
    else {
      break;
    }
  }
}

void CholeskyEliminationTree::pickRelinKeys(
  const vector<pair<Key, double>>& KeyDeltaVec,
  const std::unordered_set<Key>& allFixedKeys,
  const std::unordered_set<Key>& allMarginalizedKeys,
  int64_t remainingCycles,
  int num_threads,
  double relinThresh,
  KeySet* newRelinKeys) {

  for(sharedClique clique : cliques_) {
    if(!clique) { continue; }
    if(clique->costStatus != COST_UNMARKED || clique->nextCostStatus != COST_UNMARKED) {
      cout << "invalid cost status at clique " << *clique << " status: " << clique->costStatus << " " << clique->nextCostStatus << endl;
      exit(1);
    }
  }

  newRelinKeys->clear();

  // First get the new unmapped keys and new factors
  vector<RemappedKey> newRemappedKeys;
  for(int index = lastNumKeys - 1; index < orderingToKey_.size() - 1; index++) {
    newRemappedKeys.push_back(orderingToKey_[index]);
  }
  lastNumKeys = orderingToKey_.size();
  vector<sharedFactorWrapper> newFactors;
  for(int index = lastNumFactors; index < factors_.size(); index++) {
    newFactors.push_back(factors_[index]);
  }
  lastNumFactors = factors_.size();

  
  struct IsNewKey {
    const vector<RemappedKey>& newRemappedKeys;
    IsNewKey(const vector<RemappedKey>& newRemappedKeys_in) 
      : newRemappedKeys(newRemappedKeys_in) {}

    bool operator()(const RemappedKey& remappedKey) const {
      for(const Key& newRemappedKey : newRemappedKeys) {
        if(newRemappedKey == remappedKey) { return true; }
      }
      return false;
    }

  };

  struct IsFixedKey {
    const unordered_set<Key>& allFixedKeys;
    IsFixedKey(const unordered_set<Key>& allFixedKeys_in) 
      : allFixedKeys(allFixedKeys_in) {}

    bool operator()(const Key& key) const {
      return allFixedKeys.find(key) != allFixedKeys.end();
    }
  };

  struct IsMarginalizedKey {
    const unordered_set<Key>& allMarginalizedKeys;
    IsMarginalizedKey(const unordered_set<Key>& allMarginalizedKeys_in) 
      : allMarginalizedKeys(allMarginalizedKeys_in) {}

    bool operator()(const Key& key) const {
      return allMarginalizedKeys.find(key) != allMarginalizedKeys.end();
    }
  };

  IsNewKey isNewKey(newRemappedKeys);
  IsFixedKey isFixedKey(allFixedKeys);
  IsMarginalizedKey isMarginalizedKey(allMarginalizedKeys);

  // 1. First compute the cost of new factors
  //    This is essentially the cost of the AtA of the new factors
  //    and the of the cholesky from the connected keys to the root
  set<RemappedKey> affectedOldRemappedKeys;
  for(sharedFactorWrapper factor : newFactors) {
    for(RemappedKey remappedKey : factor->remappedKeys()) {
      if(remappedKey != 0 && !isNewKey(remappedKey)) {
        affectedOldRemappedKeys.insert(remappedKey);
      }
    }
  }

  checkCliqueParallelizable();

  vector<sharedClique> allUpdatedCliques;

  int64_t new_factor_cost = 0;

  for(RemappedKey remappedKey : affectedOldRemappedKeys) {
    vector<sharedClique> updatedCliques;
    int64_t cost = computeCost(remappedKey, num_threads, &updatedCliques);
    new_factor_cost += cost;

    commitCost(updatedCliques, &allUpdatedCliques);
  }

  // Compute cost of backsolve
  int64_t backsolve_cost = 0;
  for(sharedNode node : nodes_) {
    if(!node) { continue; }

    if(isNewKey(node->key)) { continue; }

    Key unmappedKey = unmapKey(node->key);

    if(isMarginalizedKey(unmappedKey)) { continue; }

    sharedClique clique = node->clique();
    if(clique->frontKey() != node->key) { continue; }
    if(node->key == 0 || clique == root_) { continue; }
    backsolve_cost += clique->computeCostBacksolve(num_threads);
  }

  remainingCycles -= new_factor_cost + backsolve_cost;

  cout << "backsolve cost = " << backsolve_cost << " new_factor_cost = " << new_factor_cost << " remainingCycles = " << remainingCycles << endl;

  if(remainingCycles <= 0) { 
    cout << "no relin possible!" << endl;
    resetCost(allUpdatedCliques);
    return; 
  }


  vector<pair<Key, double>> newKeyDeltaVec;
  newKeyDeltaVec.reserve(KeyDeltaVec.size());

  // static double force_thresh = 3.5;
  static double force_thresh = 100000000;
  const double min_force_thresh = 2;
  const double max_force_thresh = 3.5;

  int num_relin_keys = 0;
  int num_min_force_thresh_keys = 0;
  int min_force_keys_thresh = 10;

  for(auto& p : KeyDeltaVec) {
    Key key = p.first;
    RemappedKey remappedKey = getRemapKey(key);

    if(isNewKey(remappedKey)) { continue; }
    if(isFixedKey(key)) { continue; }
    if(isMarginalizedKey(key)) { continue; }

    if(p.second > relinThresh) {
      newKeyDeltaVec.push_back(p);
      num_relin_keys++;
    }
    if(p.second > min_force_thresh) {
      num_min_force_thresh_keys++;
    }
  }

  // double force_rate = num_relin_keys == 0? 0 : num_min_force_thresh_keys / (double) num_relin_keys;
  // if(force_rate > 0.1) {
  if(num_min_force_thresh_keys > 10) {
    force_thresh = max(force_thresh - 0.5, min_force_thresh);
  }
  else {
    force_thresh = min(force_thresh + 0.5, max_force_thresh);
  }

  int key_count = 0;
  int relin_count = KeyDeltaVec.size();

  // Sort the relevance in reverse order
  sort(newKeyDeltaVec.begin(), newKeyDeltaVec.end(), 
      [](const auto& lhs, const auto& rhs) { return lhs.second > rhs.second; });

  // 1. Go down the list of keys
  //    For each key,
  //      Compute the cost of each key
  vector<CostStatus> curRelinCostStatus(nodes_.size(), COST_UNMARKED);
  double highestUnpickedDelta = 0;
  int64_t total_relin_cost = 0;
  for(auto& p : newKeyDeltaVec) {
    Key& key = p.first;
    RemappedKey remappedKey = getRemapKey(key);

    sharedNode node = nodes_[remappedKey];

    int64_t relin_cost = computeCostRelin(remappedKey, num_threads, curRelinCostStatus);
    // cout << "relin_cost = " << relin_cost << endl;
    int64_t cost = relin_cost;
    vector<sharedClique> updatedCliques;
    for(sharedFactorWrapper factorWrapper : node->factors) {
      RemappedKey lowestKey = factorWrapper->lowestKey();
      cost += computeCost(lowestKey, num_threads, &updatedCliques);

    }

    if(cost <= remainingCycles || p.second >= force_thresh) {
      remainingCycles -= cost;
      newRelinKeys->insert(key);

      curRelinCostStatus[remappedKey] = COST_RELIN;
      commitCost(updatedCliques, &allUpdatedCliques);
      total_relin_cost += relin_cost;

    }
    else {
      uncommitCost(updatedCliques, &allUpdatedCliques);

      if(highestUnpickedDelta == 0) {
        highestUnpickedDelta = p.second;
      }
    }
    // cout << "remainingCycles = " << remainingCycles << endl;
  }

  int64_t totalCost = 0;
  for(RemappedKey key : orderingToKey_) {
    sharedClique clique = cliques_[key];
    if(!clique) { continue; }

    Key unmappedKey = unmapKey(key);

    if(clique->frontKey() != key) { continue; }
    if(isMarginalizedKey(unmappedKey)) { continue; }

    if(clique->costStatus == COST_MARKED) {
      totalCost += clique->markedCost;
    }
    else if(clique->costStatus == COST_FIXED) {
      totalCost += clique->fixedCost;
    }
    totalCost += clique->backsolveCost;
  }

  resetCost(allUpdatedCliques);
  resetCliqueParallelizable();

  for(sharedClique clique : cliques_) {
    if(!clique) { continue; }
    if(clique->costStatus != COST_UNMARKED || clique->nextCostStatus != COST_UNMARKED) {
      cout << "invalid cost status at clique " << *clique << " status: " << clique->costStatus << " " << clique->nextCostStatus << endl;
      exit(1);
    }
  }

  cout << "Remaining cycles: " << remainingCycles 
       << " Total relin keys: " << newKeyDeltaVec.size() 
       << " num relin: " << newRelinKeys->size() 
       << " relin_cost " << total_relin_cost 
       <<  " highest unpicked delta: " << highestUnpickedDelta 
       << " force_thresh: " << force_thresh 
       << " num_min_force_thresh_keys: " << num_min_force_thresh_keys 
       << " total_cost: " << totalCost << endl;

  cout << "Relin keys: ";
  for(auto k : *newRelinKeys) {
    cout << k << " ";
  }
  cout << endl;

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

  datasetgen_relinKeys.clear();
  for(const Key unmappedRelinKey : relinKeys) {
    RemappedKey relinKey = getRemapKey(unmappedRelinKey);
    datasetgen_relinKeys.push_back(relinKey);
  }

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
                                          factorIndex, newFactorIndex, factor, nullptr, this);

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
      cout << "factor has marginalized keys" << endl;
      exit(1);
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

  // Set this for staleness check in LRU
  for(const Key unmappedRelinKey : relinKeys) {
    RemappedKey remappedKey = getRemapKey(unmappedRelinKey);
    nodes_[remappedKey]->lastRelinStep = keyTransformMap_.size();
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
  // datasetgen: when marking a key, mark all the keys in the clique
  RemappedKey curKey = curClique->frontKey();
  // RemappedKey curKey = key;
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

  // if(postOrder_) {
  // In the branched relaxed supernode case, must post order every time
  if(true) {
    // sortedMarkedKeys in this case will be exactly the same as partialOrdering
    postReordering(sortedMarkedKeys);
  }

  allocateStack();

// #ifdef DEBUG
  checkInvariant_afterSymbolic();
// #endif
}

void CholeskyEliminationTree::symbolicEliminateKey(const RemappedKey key) {

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

  clique->populateBlockIndices(colStructure);

  // if(clique->nodes.front()->key != 0 &&
  //     clique->children.size() == 1) {
  //   // Clique 0 cannot be merged and will always be root
  //   sharedClique childClique = *(clique->children.begin());
  //   if(childClique->marked() 
  //       && colStructure.size() == childClique->blockIndices.size() - childClique->cliqueSize()) {
  //     clique->mergeClique2(childClique);
  //     // Need to update clique pointer to current clique, which is the child clique
  //     assert(childClique.unique());
  //     childClique.reset();

  //   }
  //   // gemmini-integrate: trying to merge more cliques together
  //   else if(childClique->marked() 
  //       && childClique->cliqueSize() < 4) {
  //     clique->mergeClique2(childClique);
  //     // Need to update clique pointer to current clique, which is the child clique
  //     assert(childClique.unique());
  //     childClique.reset();
  //   }
  // }
  
  if(clique->nodes.front()->key != 0) {
    // Clique 0 cannot be merged and will always be root
    // datasetgen: Try to merge more cliques together into relaxed supernodes
    vector<sharedClique> childCliques;
    childCliques.insert(childCliques.begin(), 
                        clique->children.begin(),
                        clique->children.end());
    for(sharedClique childClique : childCliques) {
      if(!childClique->marked()) { continue; }

      bool mergeFlag = false;
      const int relaxed_supernode_min_size = 24;
      if(clique->blockIndices.size() == childClique->blockIndices.size() - childClique->cliqueSize()) { 
        mergeFlag = true; 
      }
      else if(get<BLOCK_INDEX_ROW>(childClique->blockIndices[childClique->cliqueSize()]) < relaxed_supernode_min_size)  {
        mergeFlag = true;

      }
      else {
      }
      if(mergeFlag) {
        clique->mergeClique2(childClique);
        childClique.reset();
      }
    }
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
  // cout << "Ordering: ";
  // for(auto k : *partialOrdering) {
  //   cout << k << " ";
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

  for(RemappedKey k : partialOrdering) {
    assert(cliques_[k] != nullptr);
  }

  // Adjust ordering_. Shift not reordered keys to the front
  auto minIt = min_element(partialOrdering.begin(), partialOrdering.end(), orderingLess_);
  RemappedKey minReorderedKey = *minIt;
  size_t minReorderedIndex = keyToOrdering_[minReorderedKey];
  size_t fixedIndex = minReorderedIndex;
  vector<RemappedKey> newOrderingToKey(orderingToKey_);
  for(size_t i = minReorderedIndex; i < orderingToKey_.size(); i++) {
    RemappedKey key = orderingToKey_[i];
    sharedClique clique = cliques_[key];

    if(!clique || !cliques_[key]->marked()) {
      newOrderingToKey[fixedIndex] = key;
      fixedIndex++;
    }
  }
  
  // Add newly reordered keys to the back of the ordering
  if(fixedIndex + partialOrdering.size() != orderingToKey_.size()) {
    cout << fixedIndex << " " << partialOrdering.size() << " " << orderingToKey_.size() << endl;
  }
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

  assert(orderingToKey_.back() == 0);
  assert(keyToOrdering_[0] == keyToOrdering_.size() - 1);

  for(int i = 0; i < orderingToKey_.size(); i++) {
    auto key = orderingToKey_[i];
    if(keyToOrdering_[key] != i) {
      cout << key << " " << keyToOrdering_[key] << " " << i << endl;
    }
    assert(keyToOrdering_[key] == i);
  }

}

void CholeskyEliminationTree::postReordering(const vector<RemappedKey>& partialOrdering) {
  vector<RemappedKey> newOrderingToKey(orderingToKey_);

  RemappedKey minReorderedKey = orderingToKey_[0]; // partialOrdering.front();
  int minReorderedIndex = keyToOrdering_[minReorderedKey];

  assert(keyToOrdering_[0] == keyToOrdering_.size() - 1);
  assert(orderingToKey_.back() == 0);
  assert(partialOrdering.back() == 0);

  // If come across the last key of a clique, reorder all keys in the clique to be 
  // right in front of it. This is needed for branched relaxed supernodes
  int curIndex = orderingToKey_.size() - 1; 
  int removedKeyIndex = 0;
  for(int i = curIndex; i >= minReorderedIndex; i--) {
    RemappedKey key = orderingToKey_[i];
    sharedClique clique = cliques_[key];
    
    if(!clique) {
      newOrderingToKey[removedKeyIndex] = key;
      removedKeyIndex++;
      continue;
    }
    else if(key == clique->backKey()) {
      for(int k = clique->cliqueSize() - 1; k >= 0; k--) {
        RemappedKey cliqueKey = get<BLOCK_INDEX_KEY>(clique->blockIndices[k]);
        assert(curIndex < newOrderingToKey.size());
        newOrderingToKey[curIndex] = cliqueKey;
        curIndex--;
      }
    }
  }

  assert(removedKeyIndex == curIndex + 1);

  orderingToKey_ = newOrderingToKey;

  for(size_t i = minReorderedIndex; i < orderingToKey_.size(); i++) {
    RemappedKey key = orderingToKey_[i];
    keyToOrdering_[key] = i;
  }

  for(int i = 0; i < orderingToKey_.size(); i++) {
    auto key = orderingToKey_[i];
    if(keyToOrdering_[key] != i) {
      cout << key << " " << keyToOrdering_[key] << " " << i << endl;
    }
    assert(keyToOrdering_[key] == i);
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
          // nodes_[key]->clique()->setStatusEdit();
          nodes_[key]->clique()->setStatusReconstruct();
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
    // Gather columns from wherever they are living now
    // Need to do a scatter operation based on the old blockIndices on to the new one
    // The clique who owns the data shall be responsible for the blockIndices
    GEMMINI_TYPE* matrixData = workspace_.get_ptrs(clique->workspaceIndex)[0];
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

  auto B = Eigen::Block<Eigen::Map<ColMajorMatrix<GEMMINI_TYPE>>>
            (m, diagWidth, 0, subdiagHeight, diagWidth);
  auto C = Eigen::Block<Eigen::Map<ColMajorMatrix<GEMMINI_TYPE>>>
            (m, diagWidth, diagWidth, subdiagHeight, subdiagHeight);

  assert(get<BLOCK_INDEX_KEY>(blockIndices.back()) == 0);

  bool processGrouped = true;
  size_t destSize = destCols.size();
  size_t totalSize = blockIndices.size() - 1;
  // if(destSize * 2 < totalSize) {
  //     processGrouped = false;
  // }

  if(processGrouped) {
    // We want to get C^T = B^T^T B^T, we have column major B, which is row major B^T
    syrk(subdiagHeight, subdiagHeight, diagWidth,
           &m(diagWidth, 0), &m(diagWidth, 0), &m(diagWidth, diagWidth),
           totalHeight, totalHeight, totalHeight,
           sign, 1, 
           true, false);

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

  Eigen::Map<ColMajorMatrix<GEMMINI_TYPE>> m = workspace_.get_matrices(clique->workspaceIndex)[0];

  size_t totalHeight = clique->height();
  size_t diagWidth = clique->width();

  auto DB = Eigen::Block<Eigen::Map<ColMajorMatrix<GEMMINI_TYPE>>>(m, 0, 0, totalHeight, diagWidth);
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
struct CheckEditReconstructNew {
  bool operator()(const MarkedStatus& status) const { 
    return (status == EDIT) || (status == RECONSTRUCT) || (status == NEW);
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
    assert(nodes_[key]->status() != FIXED);
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
        // assert(factorWrapper->getLinearFactorType() != NONE);

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
        factorWrapper->updateHessian(m, 1, info, CheckEditReconstructNew());
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
  // cout << "block indices: ";
  // clique->printBlockIndices(cout);
  // cout << endl;
  // m.triangularView<Eigen::StrictlyUpper>().setZero();
  // // cout << "Before eliminate. m = \n" << block(m, 0, 0, totalHeight, bWidth) << endl << endl;

  // // cout << "Before eliminate" << endl;
  // // for(int j = 0; j < bWidth; j += 3) {
  // //     cout << "j = " << j << endl;
  // //     cout << "m = \n" << block(m, 0, j, totalHeight, 3) << endl << endl;
  // // }

  Eigen::Block<Eigen::Map<ColMajorMatrix<GEMMINI_TYPE>>> D = block(m, 0, 0, bWidth, bWidth);
  Eigen::LLT<Eigen::Ref<Eigen::Block<Eigen::Map<ColMajorMatrix<GEMMINI_TYPE>>>>> llt(D);
  if(llt.info() == Eigen::NumericalIssue) {
    cout << "Diagonal block not positive definite!" << endl;
    clique->printClique(cout);
    cout << "D = \n" << D << endl << endl;
    exit(1);
  }

  auto L = D.triangularView<Eigen::Lower>();
  auto B = block(m, bWidth, 0, bHeight, bWidth);
  auto C = block(m, bWidth, bWidth, bHeight, bHeight);

  L.solveInPlace(B.transpose());

  if(bHeight != -1) {
    // We want to get C^T = B^T^T B^T, we have column major B, which is row major B^T
    syrk(bHeight, bHeight, bWidth,
           &m(bWidth, 0), &m(bWidth, 0), &m(bWidth, bWidth),
           totalHeight, totalHeight, totalHeight,
           -1, 1, 
           true, false);
  }

  // // DEBUG
  // cout << "Clique: " << *clique << endl << endl;
  // cout << "block indices: ";
  // clique->printBlockIndices(cout);
  // cout << endl;
  // m.triangularView<Eigen::StrictlyUpper>().setZero();
  // // cout << "After eliminate. m = \n" << block(m, 0, 0, totalHeight, bWidth) << endl << endl;
  // // cout << "After eliminate" << endl;
  // // for(int j = 0; j < bWidth; j += 3) {
  // //     cout << "j = " << j << endl;
  // //     cout << "m = \n" << block(m, 0, j, totalHeight, 3) << endl << endl;
  // // }
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
    std::shared_ptr<vector<GEMMINI_TYPE>> matrixSource 
      = std::make_shared<vector<GEMMINI_TYPE>>(r * c);
    std::shared_ptr<BlockIndexVector> blockIndicesSource 
      = std::make_shared<BlockIndexVector>(clique->blockIndices);

    // // DEBUG
    // auto m = workspace_.get_matrices(clique->workspaceIndex)[0];
    // cout << "m rows = " << m.rows() << " m cols = " << m.cols() << " r = " << r << " c = " << c << endl;
    // assert(r == m.rows() && c <= m.cols());

    memcpy(matrixSource->data(), m_ptr, sizeof(GEMMINI_TYPE) * r * c);

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

// #ifdef DEBUG
  checkInvariant_afterBackSolve();
// #endif
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
  GemminiVector delta = block(m, totalHeight - 1, 0, 1, diagWidth).transpose();

  if(subdiagHeight > 1) {
    Eigen::Matrix<GEMMINI_TYPE, Eigen::Dynamic, 1> gatherX(subdiagHeight - 1);

    // Gather deltas this clique depends on, don't need last row
    for(int i = cliqueSize; i < blockIndices.size() - 1; i++) {
      auto[key, row, height] = blockIndices[i];
      Key unmappedKey = unmapKey(key);

      row -= diagWidth;

      gatherX.block(row, 0, height, 1) = delta_ptr->at(unmappedKey).cast<GEMMINI_TYPE>();
    }

    gemv(diagWidth, subdiagHeight - 1, 
         &m(diagWidth, 0), &gatherX(0), &delta(0), 
         totalHeight, 
         -1);

    // auto B = block(m, diagWidth, 0, subdiagHeight - 1, diagWidth); // sub-diagonal blocks
    // delta -= B.transpose() * gatherX;

    // cout << "delta = " << delta << endl;;
    // cout << "delta copy = " << delta_copy << endl;
  }
  
  // Solve diagonal
  auto D = block(m, 0, 0, diagWidth, diagWidth);
  auto LT = D.transpose().triangularView<Eigen::Upper>();
  LT.solveInPlace(delta);

  Vector delta_double = delta.cast<double>();

  // Copy delta back into delta_ptr
  for(int i = 0; i < cliqueSize; i++) {
    auto[key, row, height] = blockIndices[i];
    Key unmappedKey = unmapKey(key);

    Vector delta_diff = delta_double.block(row, 0, height, 1) - delta_ptr->at(unmappedKey);
    delta_ptr->at(unmappedKey) = delta_double.block(row, 0, height, 1);

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

void CholeskyEliminationTree::selectStaleSubtree(
        const unordered_set<Key>& allFixedKeys,
        const unordered_set<Key>& allMarginalizedKeys,
        int targetTreeSize, 
        vector<Key>* marginalizedKeys,
        vector<Key>* fixedKeys) {

  marginalizedKeys->clear();
  fixedKeys->clear();

  vector<pair<sharedClique, bool>> stack(1, {root_, false});

  size_t total_size = 0;

  // First compute the subtree sizes of each clique
  while(!stack.empty()) {
    auto& curPair = stack.back();
    sharedClique curClique = curPair.first;
    bool& expanded = curPair.second;

    if(!expanded) {
      expanded = true;    

      for(sharedClique childClique : curClique->children) {
        if(allMarginalizedKeys.find(childClique->backKey()) != allMarginalizedKeys.end()) {
          continue;
        }
        stack.push_back({childClique, false});
      }
    }
    else {
      stack.pop_back();

      curClique->memSize = curClique->width() * curClique->width();
      curClique->subdiagSize = curClique->subdiagonalHeight() * curClique->width();

      // for(sharedNode node : curClique->nodes) {
      //   for(auto factor : node->factors) {
      //     if(node->key != factor->lowestKey()) { continue; }
      //     int total_width = 0;
      //     for(const auto&[key, col, width] : factor->blockIndices()) {
      //       total_width += width;
      //     }
      //     const auto& cachedMatrix = factor->getCachedMatrix();
      //     curClique->memSize += cachedMatrix.rows() * total_width;
      //   }
      // }

      for(sharedClique childClique : curClique->children) {
        if(allMarginalizedKeys.find(childClique->backKey()) != allMarginalizedKeys.end()) {
          continue;
        }
        curClique->memSize += childClique->memSize + childClique->subdiagSize;
      }

      curClique->lastRelinStep = 0;
      for(sharedNode node : curClique->nodes) {
        if(node->lastRelinStep > curClique->lastRelinStep) {
          curClique->lastRelinStep = node->lastRelinStep;
        }
      }

      for(sharedClique childClique : curClique->children) {
        if(childClique->lastRelinStep > curClique->lastRelinStep) {
          curClique->lastRelinStep = childClique->lastRelinStep;
        }
      }
    }
  }

  int totalMemSize = root_->memSize;

  cout << "totalMemSize = " << totalMemSize << " targetTreeSize = " << targetTreeSize << endl;

  if(totalMemSize <= targetTreeSize) {
    return;
  }

  int discardSize = totalMemSize - targetTreeSize;
  
  cout << "discardSize = " << discardSize << endl;

  stack.push_back({root_, false});

  vector<sharedClique> candidates;

  // First compute the subtree sizes of each clique
  while(!stack.empty()) {
    auto& curPair = stack.back();
    stack.pop_back();
    sharedClique curClique = curPair.first;
    bool& expanded = curPair.second;

    // curClique is assumed to have memSize larger than discardSize
    // If child has memSize larger than discardSize, expand on that child
    // If no child or no child is large than discardSize, add curClique to candidates

    if(curClique->children.empty()) {
        candidates.push_back(curClique);
    }
    else {
      bool flag = false;
      for(sharedClique childClique : curClique->children) {
        if(childClique->memSize >= discardSize) {
          flag = true;
          stack.push_back({childClique, false});
        }
      }
      if(!flag) {
        if(curClique == root_) {
          cout << "Discarding root!" << endl;
          cout << "First child clique" << *(*curClique->children.begin()) << (*curClique->children.begin())->memSize << endl;
          exit(1);
        }
        candidates.push_back(curClique);
      }
    }
  }

  struct CandidateLess {
    bool operator()(const sharedClique& lhs, const sharedClique& rhs) const {
      return lhs->lastRelinStep == rhs->lastRelinStep? lhs->memSize < rhs->memSize :
              lhs->lastRelinStep < rhs->lastRelinStep;
    }
  };

  sort(candidates.begin(), candidates.end(), CandidateLess());

  sharedClique selClique = candidates.front();

  stack.push_back({selClique, false});

  while(!stack.empty()) {
    auto& curPair = stack.back();
    stack.pop_back();
    sharedClique curClique = curPair.first;
    bool& expanded = curPair.second;

    for(int i = 0; i < curClique->cliqueSize(); i++) {
      RemappedKey remappedKey = get<BLOCK_INDEX_KEY>(curClique->blockIndices[i]);
      marginalizedKeys->push_back(unmapKey(remappedKey));
    }

    for(sharedClique childClique : curClique->children) {
      stack.push_back({childClique, false});
    }
  }

  for(int i = selClique->cliqueSize(); i < selClique->blockIndices.size() - 1; i++) {
    RemappedKey remappedKey = get<BLOCK_INDEX_KEY>(selClique->blockIndices[i]);
    fixedKeys->push_back(unmapKey(remappedKey));
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

  cout << "marginalize. should not be here" << endl;

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
  factorMatrix.selfadjointView().rankUpdate(B.cast<double>(), -1);

  // Make HessianFactor
  HessianFactor::shared_ptr marginalFactor = boost::make_shared<HessianFactor>(unmappedKeys, factorMatrix);

  size_t factorIndex = factors_.size();
  sharedFactorWrapper factorWrapper = std::make_shared<FactorWrapper>(
                                        factorIndex, -1, nullptr, marginalFactor, this);

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

// Reimplement marginalization. We should just remove all traces of the marginalized keys
// Including factors
void CholeskyEliminationTree::marginalizeLeaves2(
    const std::vector<Key>& leafKeys,
    boost::optional<FactorIndices&> marginalFactorsIndices,
    boost::optional<FactorIndices&> deletedFactorsIndices) {

  // First unmap the leaf keys and reorder
  vector<RemappedKey> remappedLeafKeys;
  for(Key unmappedKey : leafKeys) {
    RemappedKey key = getRemapKey(unmappedKey);
    remappedLeafKeys.push_back(key);
  }

  sort(remappedLeafKeys.begin(), remappedLeafKeys.end(), orderingLess_);

  // Mark all the cliques and nodes that will be marginalized
  vector<sharedClique> marginalizedCliques;
  for(RemappedKey key : remappedLeafKeys) {
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

  vector<sharedFactorWrapper> deletedFactors;
  vector<sharedFactorWrapper> marginalFactors;

  for(sharedClique clique : marginalizedCliques) {
    marginalizeClique2(clique, &deletedFactors, &marginalFactors);
  }

  // Remove all traces of marginalized keys. Delete cliques, nodes, keys, and factors
  for(RemappedKey key : remappedLeafKeys) {
    nodes_[key]->deleteNode();
    nodes_[key] = nullptr;
    unmappedKeys_[key] = -1;
    // keyToOrdering_[key] = -1;
  }

  for(sharedClique clique : marginalizedCliques) {
    cout << "deleting clique: " << *clique << endl;
    clique->deleteClique();
  }

  for(sharedFactorWrapper factor : deletedFactors) {
    // TODO: Add to deletedFactorsIndices
    // Not killing this because we are printing out factors in extract key
    for(RemappedKey key : factor->remappedKeys()) {
      if(key == 0) { continue; }
      if(nodes_[key]) {
        cout << "node " << key << " is supposed to be null!" << endl;
        exit(1);
      }
    }
    factor->setStatusRemoving();
    factors_[factor->factorIndex()] = nullptr;
  }

  for(sharedFactorWrapper factor : factors_) {
    if(!factor) { continue; }
    for(RemappedKey key : factor->remappedKeys()) {
      if(key == 0) { continue; }
      if(!nodes_[key]) {
        cout << "Factor " << *factor << endl;
        cout << "node " << key << " is not supposed to be null!" << endl;
        exit(1);
      }
    }
  }

  // Don't do anything to orderingToKey_ and keyToOrdering_. Too complicated

  int64_t total_size = 0;
  vector<sharedClique> stack(1, root_);
  while(!stack.empty()) {
    sharedClique curClique = stack.back();
    stack.pop_back();

    total_size += curClique->width() * curClique->height();

    // for(sharedNode node : curClique->nodes) {
    //   for(auto factor : node->factors) {
    //     if(node->key != factor->lowestKey()) { continue; }
    //     int total_width = 0;
    //     for(const auto&[key, col, width] : factor->blockIndices()) {
    //       total_width += width;
    //     }
    //     const auto& cachedMatrix = factor->getCachedMatrix();
    //     total_size += cachedMatrix.rows() * total_width;
    //   }
    // }

    for(sharedClique childClique : curClique->children) {
      stack.push_back(childClique);
    }
  }

  cout << "totalMemSize after marginalization: " << total_size << endl;

#ifdef DEBUG
  checkInvariant_afterMarginalize();
#endif
}

void CholeskyEliminationTree::marginalizeClique2(
      sharedClique clique,
      vector<sharedFactorWrapper>* deletedFactors,
      vector<sharedFactorWrapper>* marginalFactors) {

  // For now, for LRU, assume that we cannot split cliques
  if(nodes_[clique->backKey()]->status() != MARGINALIZED) {
    cout << "Marginalization within cliques not implemented!" << endl;
    exit(1);
  }

  // All of this clique's children should be marginalized
  // this is required of the input of leafKeys
  for(sharedClique childClique : clique->children) {
    if(childClique->status() != MARGINALIZED) {
      cout << "Marginalizing clique that has unmarginalized children" << endl;
      exit(1);
    }
  }

  cout << "Marginalizing clique " << *clique << endl;
  
  // Set all factors in each node to be LINEAR
  // and remove all keys and data that are marginalized out
  // If all the keys in a factor is marginalized out, remove that factor
  for(sharedNode node : clique->nodes) {
    if(node->status() != MARGINALIZED) { continue; }
    for(sharedFactorWrapper factor : node->factors) {
      if(node->key == factor->lowestKey()) {  // If node is not lowest key, then factor would have been deleted by someone else
        bool deleteFactor = factor->marginalizeKeys();
        if(deleteFactor) {
          deletedFactors->push_back(factor);
        }
        else {
          assert(factor->status() == LINEAR);
        }
      }
    }
  }

  // If a clique has not connections to an unmarginalized node
  // we can just return after deleting factors
  RemappedKey secondToLastKey 
    = get<BLOCK_INDEX_KEY>(clique->blockIndices[clique->blockIndices.size() - 2]);
  if(nodes_[secondToLastKey]->status() == MARGINALIZED) {
    return;
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
  factorMatrix.selfadjointView().rankUpdate(B.cast<double>(), -1);

  // Make HessianFactor
  HessianFactor::shared_ptr marginalFactor = boost::make_shared<HessianFactor>(unmappedKeys, factorMatrix);

  size_t factorIndex = factors_.size();
  sharedFactorWrapper factorWrapper = std::make_shared<FactorWrapper>(
                                        factorIndex, -1, nullptr, marginalFactor, this);

  factors_.push_back(factorWrapper);

  for(RemappedKey key : factorWrapper->remappedKeys()) {
    if(key == 0) { continue; }
    nodes_[key]->addFactor(factorWrapper);
    // nodes_[key]->setStatusFixed();
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

void CholeskyEliminationTree::extractSubtree(std::ostream& os, int size) const {
  queue<sharedClique> q;
  vector<sharedClique> collectedCliques;
  set<sharedFactorWrapper> collectedFactors;
  q.push(root_);

  int collectedSize = 0;
  int numNodes = 0;

  while(!q.empty()) {
    sharedClique clique = q.front();
    q.pop();

    collectedCliques.push_back(clique);

    for(sharedNode node : clique->nodes) {
      collectedSize += node->width;
      numNodes++;

      for(sharedFactorWrapper factorWrapper : node->factors) {
        if(factorWrapper->lowestKey() == node->key) {
          collectedFactors.insert(factorWrapper);
        }
      }
    }

    if(collectedSize >= size) {
      break;
    }

    for(sharedClique child : clique->children) {
      q.push(child);
    }
  }

  os << "factors\n" << collectedFactors.size() << endl;
  for(sharedFactorWrapper factorWrapper : collectedFactors) {
    cout << *factorWrapper << endl;
    os << factorWrapper->factorIndex() << " " << factorWrapper->getCachedMatrix().rows() << " " << factorWrapper->getCachedMatrix().cols() << " " << factorWrapper->remappedKeys().size() << " ";
    for(RemappedKey remappedKey : factorWrapper->remappedKeys()) {
      os << remappedKey << " ";
    }
    os << endl;
  }

  unordered_map<sharedClique, int> cliqueMap;

  int c = 0;
  for(auto it = collectedCliques.rbegin(); it != collectedCliques.rend(); it++, c++) {
    cout << *it << " " << c << endl;
    cliqueMap.insert({*it, c});
  }

  os << endl << "cliques\n" << collectedCliques.size() << endl;

  for(auto it = collectedCliques.rbegin(); it != collectedCliques.rend(); it++) {
    sharedClique clique = *it;

    os << "clique\n" << cliqueMap.at(clique) << " ";
    if(clique->parent() != nullptr) {
      os << cliqueMap.at(clique->parent());
    }
    else {
      os << "-1";
    }
    os << " ";

    os << clique->nodes.size() << " " << clique->blockIndices.size() << endl;


    for(auto& [key, row, height] : clique->blockIndices) {
      os << key << " " << row << " " << height << endl;
    }

    os << endl;

    vector<sharedFactorWrapper> cliqueFactors;

    for(sharedNode node : clique->nodes) {
      for(sharedFactorWrapper factorWrapper : node->factors) {
        if(factorWrapper->lowestKey() == node->key) {
          cliqueFactors.push_back(factorWrapper);
        }
      }
    }

    os << "clique factors\n" << cliqueFactors.size() << endl;

    for(sharedFactorWrapper factorWrapper : cliqueFactors) {
      os << factorWrapper->factorIndex() << endl;
    }

    os << endl;
  }

  os << endl;
  
}

void CholeskyEliminationTree::extractFullTree(std::ostream& os, bool print_values) const {

  os << "ordering and width" << endl;
  os << orderingToKey_.size() << endl;
  for(RemappedKey remappedKey : orderingToKey_) {
    sharedNode node = nodes_[remappedKey];
    if(!node) {
      os << remappedKey << " " << 0 << endl;
    }
    else {
      os << remappedKey << " " << colWidth(remappedKey) << endl;
    }
  }
  os << endl;

  int num_active_factors = 0;
  for(sharedFactorWrapper factorWrapper : factors_) {
    if(factorWrapper) {
      num_active_factors++;
    }
  }

  os << "factors" << endl;
  os << num_active_factors << endl;
  for(sharedFactorWrapper factorWrapper : factors_) {
    if(!factorWrapper) { continue; }
    os << factorWrapper->factorIndex() << " " 
       << factorWrapper->getCachedMatrix().rows() << " " 
       << factorWrapper->getCachedMatrix().cols() << " " 
       << factorWrapper->remappedKeys().size() << endl;
    for(RemappedKey factorKey : factorWrapper->remappedKeys()) {
      if(!nodes_[factorKey]) {
        cout << "factorKey: " << factorKey << " is NULL!" << endl;
        exit(1);
      }
      os << factorKey << " ";
    }
    os << endl;
    if(print_values) {
      os << factorWrapper->getCachedMatrix() << endl;
    }
  }
  os << endl;

  int clique_count = 0;
  for(RemappedKey remappedKey : orderingToKey_) {
    sharedNode node = nodes_[remappedKey];
    if(!node) { continue; }

    sharedClique clique = node->clique();
    if(clique->frontKey() != remappedKey) { continue; }
    clique_count++;
  }

  os << "cliques" << endl;
  os << clique_count << endl;

  for(RemappedKey remappedKey : orderingToKey_) {
    sharedNode node = nodes_[remappedKey];
    if(!node) { continue; }

    sharedClique clique = nodes_[remappedKey]->clique();
    if(clique->frontKey() != remappedKey) { continue; }

    clique_count++;

    os << clique->height() << " " << clique->width() << " "
       << clique->cliqueSize() << " " << clique->blockIndices.size() << endl;
    for(int i = 0; i < clique->blockIndices.size(); i++) {
      const auto&[key, row, height] = clique->blockIndices[i];
      os << key << " " << row << " " << height << endl;
    }

    if(clique->gatherSources.size() != 1) {
      cerr << "Clique is not well formed!" << endl;
      exit(1);
    }
    
    if(print_values) {
      auto m = clique->gatherSources.front().matrix();
      m.triangularView<Eigen::StrictlyUpper>().setZero();
      os << m << endl;
    }

  }
  os << endl;

  os << "new keys\n" << datasetgen_newKeys.size() << endl;
  for(RemappedKey k : datasetgen_newKeys) {
    os << k << endl;
  }
  os << endl;

  os << "relin keys\n" << datasetgen_relinKeys.size() << endl;
  for(RemappedKey k : datasetgen_relinKeys) {
    os << k << endl;
  }

  os << endl;


}

void CholeskyEliminationTree::extractDelta(std::ostream& os, VectorValues& delta) const {
  os << "deltas" << endl;
  os << delta.size() << endl;

  for(RemappedKey remappedKey = 0; remappedKey < unmappedKeys_.size(); remappedKey++) {
    if(remappedKey == 0) { continue; } 

    Key unmappedKey = unmapKey(remappedKey);

    os << delta.at(unmappedKey).transpose() << endl;
  }
  os << endl;
}

void CholeskyEliminationTree::extractPredictedCycles(std::ostream& os, int num_threads) {
  for(auto clique : cliques_) {
    if(!clique) { continue; }
    if(clique->costStatus != COST_UNMARKED || clique->nextCostStatus != COST_UNMARKED) {
      cout << "cost status not reset for clique: " << *clique << endl;
      exit(1);
    }
  }

  checkCliqueParallelizable();

  // Run cost model for known relin keys and new keys
  int total_predicted_cost = 0;
  vector<sharedClique> updatedCliques;
  vector<sharedClique> allUpdatedCliques;
  for(RemappedKey k : datasetgen_newKeys) {
    sharedNode node = nodes_[k];
    for(sharedFactorWrapper factor : node->factors) {
      RemappedKey lowestKey = factor->lowestKey();
      total_predicted_cost += computeCost(lowestKey, num_threads, &updatedCliques);
      commitCost(updatedCliques, &allUpdatedCliques);
    }
  }
  int relin_cost = 0;
  vector<CostStatus> curRelinCostStatus(nodes_.size(), COST_UNMARKED);
  for(RemappedKey k : datasetgen_relinKeys) {
    sharedNode node = nodes_[k];
    relin_cost += computeCostRelin(k, num_threads, curRelinCostStatus);
    curRelinCostStatus[k] = COST_RELIN;
    for(sharedFactorWrapper factor : node->factors) {
      RemappedKey lowestKey = factor->lowestKey();
      total_predicted_cost += computeCost(lowestKey, num_threads, &updatedCliques);
      commitCost(updatedCliques, &allUpdatedCliques);
    }
  }

  int clique_count = 0;
  for(RemappedKey remappedKey : orderingToKey_) {
    sharedNode node = nodes_[remappedKey];
    if(!node) { continue; }

    sharedClique clique = node->clique();
    if(clique->frontKey() != remappedKey) { continue; }
    clique_count++;

    total_predicted_cost += clique->computeCostBacksolve(num_threads);
  }

  int index = 0;

  int total = 0;

  for(RemappedKey remappedKey : orderingToKey_) {
    sharedNode node = nodes_[remappedKey];
    if(!node) { continue; }

    sharedClique clique = node->clique();
    if(clique->frontKey() != remappedKey) { continue; }

    if(clique->costStatus == COST_MARKED) {
      total += clique->symCost + clique->AtACost + clique->cholCost + clique->addCost + clique->backsolveCost;
    }
    else if(clique->costStatus == COST_FIXED) {
      total += clique->symCost + clique->AtACost + clique->syrkCost + clique->addCost + clique->backsolveCost;
    }
    else {
    }

    index++;
  }
  os << "totalCost: " << total << endl;

  resetCost(allUpdatedCliques);
  resetCliqueParallelizable();

  for(auto clique : cliques_) {
    if(!clique) { continue; }
    if(clique->costStatus != COST_UNMARKED || clique->nextCostStatus != COST_UNMARKED) {
      cout << "cost status not reset for clique: " << *clique << endl;
      exit(1);
    }
  }

  setCostReorder(false);
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
    if(!clique) { continue; }
    clique->checkInvariant();
  }

  assert(orderingToKey_.back() == 0);
  assert(keyToOrdering_[0] == keyToOrdering_.size() - 1);

  for(int i = 0; i < orderingToKey_.size(); i++) {
    auto key = orderingToKey_[i];
    if(keyToOrdering_[key] != i) {
      cout << key << " " << keyToOrdering_[key] << " " << i << endl;
    }
    assert(keyToOrdering_[key] == i);
  }

}

void CholeskyEliminationTree::checkInvariant_afterBackSolve() const {
  for(const sharedClique& clique : cliques_) {
    if(!clique) { continue; }
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
    if(!node) { continue; }
    // One in nodes_ and one in clique->nodes
    assert(node.use_count() == 2);
  }
}

void CholeskyEliminationTree::checkInvariant_afterMarginalize() const {
  cout << "check invariant after marginalized" << endl;
  {
    for(sharedNode node : nodes_) {
      if(!node) { continue; }
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
      if(!node) { continue; }
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
