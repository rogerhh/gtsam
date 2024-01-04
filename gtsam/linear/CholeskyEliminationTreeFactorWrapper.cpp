/**
* @file    CholeskyEliminationTreeClique.cpp
* @brief   A wrapper around linear/nonlinear factors in the CholeskyEliminationTree
* @author  Roger Hsiao (rogerhh)
* @date    Mar. 16, 2023
*/

#include "CholeskyEliminationTree.h"
#include "gtsam/linear/gemmini_functions.h"
#include <gtsam/linear/CholeskyEliminationTreeFactorWrapper.h>
#include <gtsam/linear/CholeskyEliminationTreeClique.h>
#include <gtsam/linear/CholeskyEliminationTreeNode.h>
#include <gtsam/linear/JacobianFactor.h>
#include <iostream>
#include <cassert>

using namespace std;

namespace gtsam {

void CholeskyEliminationTree::FactorWrapper::markAffectedKeys(RemappedKeySet* affectedKeys) {
  if(status_ == REMOVED) {
    return;
  }
  for(RemappedKey key : remappedKeys_) {
    assert(etree->nodes_[key]->status() == UNMARKED || etree->nodes_[key]->status() == NEW);
    affectedKeys->insert(key);
  }
}


bool CholeskyEliminationTree::FactorWrapper::marginalizeKeys() {
  assert(cachedLinearFactor != nullptr);

  if(etree->nodes_[highestKey_]->status() == MARGINALIZED) { 
    status_ = REMOVED;
    return true; 
  }

  setStatusLinear();

  // Just remove the block indices corresponding to the keys that are marginalized
  vector<RemappedKey> newRemappedKeys;
  BlockIndexVector newBlockIndices;
  for(auto&[key, col, width] : blockIndices_) {
    if(etree->nodes_[key]->status() != MARGINALIZED) {
      newBlockIndices.push_back({key, col, width});
      newRemappedKeys.push_back(key);
    }
  }

  blockIndices_ = newBlockIndices;
  remappedKeys_ = newRemappedKeys;

  assert(remappedKeys_.back() == 0);

  for(const RemappedKey& key : remappedKeys_) {
    assert(etree->nodes_[key]->status() == UNMARKED);
  }

  updateLowestHighestKeys();

  return false;
}

bool CholeskyEliminationTree::FactorWrapper::hasMarginalizedKeys() const {
  for(size_t i = 0; i < remappedKeys_.size() - 1; i++) {
    const RemappedKey& key = remappedKeys_[i];
    if(etree->nodes_[key]->status() == MARGINALIZED) {
      return true;
    }
  }
  return false;
}

void CholeskyEliminationTree::FactorWrapper::setCachedLinearMatrix(sharedLinearFactor linearFactor) {
  if(const JacobianFactor* jf = dynamic_cast<const JacobianFactor*>(linearFactor.get())) {
    linearFactorType_ = JACOBIAN;
    cachedLinearMatrix_ = jf->matrixObject().matrix().template cast<GEMMINI_TYPE>();
  }
  else if(const HessianFactor* hf = dynamic_cast<const HessianFactor*>(linearFactor.get())) {
<<<<<<< HEAD
    linearFactorType_ = JACOBIAN;
    cachedLinearMatrix_ = hf->info().matrix().template cast<GEMMINI_TYPE>();
=======
    linearFactorType_ = HESSIAN;
    const GemminiMatrix& hessianUpper = hf->info().matrix().template cast<GEMMINI_TYPE>();
    cachedLinearMatrix_ = GemminiMatrix(hessianUpper.selfadjointView<Eigen::Upper>());
    // cachedLinearMatrix_ = hf->info().matrix().template cast<GEMMINI_TYPE>();
  }
  else {
    std::cout << "Factor type not supported!" << std::endl;
    exit(1);
>>>>>>> fp32-v2
  }
}

// const JacobianFactor* CholeskyEliminationTree::FactorWrapper::toJacobianFactor() const {
//   return dynamic_cast<const JacobianFactor*>(cachedLinearFactor.get());
// }
// 
// const HessianFactor* CholeskyEliminationTree::FactorWrapper::toHessianFactor() const {
//   return dynamic_cast<const HessianFactor*>(cachedLinearFactor.get());
// }

ostream& operator<<(ostream& os, const CholeskyEliminationTree::FactorWrapper& factorWrapper) {
  os << "Factor: " << factorWrapper.factorIndex() 
     << ", status: " << factorWrapper.status() 
     << ", keys: ";
  factorWrapper.printKeys(os);
  return os;
}

void CholeskyEliminationTree::FactorWrapper::printKeys(std::ostream& os) const {
  for(auto key : remappedKeys_) {
    os << key << " ";
  }
}

} // namespace gtsam
