/**
* @file    CholeskyEliminationTreeFactorWrapper.h
* @brief   A wrapper class around a Factor in the CholeskyEliminationTree. The Factor can be linear or nonlinear. The main purpose for this class is to provide a convenient lookup for remapped keys and to keep track of the state of each factor
* @author  Roger Hsiao (rogerhh)
* @date    Jun. 20, 2023
*/

#pragma once

#include <gtsam/linear/CholeskyEliminationTree.h>
#include <gtsam/linear/JacobianFactor.h>
#include <gtsam/nonlinear/NonlinearFactor.h>
#include <utility>

#include <gtsam/linear/gemmini_functions.h>

namespace gtsam {

class CholeskyEliminationTree::FactorWrapper :
  public std::enable_shared_from_this<CholeskyEliminationTree::FactorWrapper> {

public:

  sharedFactorWrapper get_ptr() {
    return shared_from_this();
  }

  typedef Eigen::Block<const Matrix> constBlock;

  CholeskyEliminationTree* etree = nullptr;
  size_t factorIndex_;
  size_t gtsamFactorIndex_;

  // blockIndices_ store the list of (key, col, width) of the block row
  // This DOES NOT change after reordering keys. Only needs updating after
  // adding new keys
  BlockIndexVector blockIndices_;

  // Dimension is not supported as the factor can be a Hessian factor 
  // and cannot easily get dim
  // size_t factorDim;

  // State of the factor. 
  // UNLINEARIZED is the state of all new nonlinear factors
  // RELINEARIZE is the state of old factors that need relienarization, we need a separate state for this in the case of edits
  // LINEARIZED is the state of nonlinear factors after linearization
  // LINEAR is the state of linear/marginalized factors
  // REMOVING factors should be subtracted from the info matrix if the clique is EDIT, but 
  // should not be added back in
  // REMOVED factors is the final state of REMOVING factors. We keep the factors for
  // now to retain flexibility
  // enum Status {UNLINEARIZED, LINEARIZED, RELINEARIZE, LINEAR, REMOVING, REMOVED};

  FactorStatus status_ = UNLINEARIZED;

  // A nonlinear factor linearizes to a linear augmented block row matrix
  // If a factor is already linear, then we instantiate directly as an augmented matrix
  sharedFactor nonlinearFactor_ = nullptr;

  // ColMajorMatrix cachedLinearFactor;
  // sharedLinearFactor cachedLinearFactor = nullptr;
  enum LinearFactorType { NONE, JACOBIAN, HESSIAN }; 
  LinearFactorType linearFactorType_ = NONE;
  ColMajorMatrix<GEMMINI_TYPE> cachedLinearMatrix_;

  RemappedKeyVector remappedKeys_;

  RemappedKey lowestKey_, highestKey_;

  // Construct a wrapper around a factor
  // If nonlinearFactor_in == nullptr, then we are constructing a linear factor
  FactorWrapper(size_t factorIndex_in,
                FactorIndex gtsamFactorIndex_in,
                sharedFactor nonlinearFactor_in,
                sharedLinearFactor cachedLinearFactor_in,
                CholeskyEliminationTree* etree_in) 
  : etree(etree_in),
    factorIndex_(factorIndex_in),
    gtsamFactorIndex_(gtsamFactorIndex_in),
    nonlinearFactor_(nonlinearFactor_in) { //,
    // cachedLinearFactor(cachedLinearFactor_in) {
    if(nonlinearFactor_ == nullptr) {

      assert(cachedLinearFactor_in != nullptr);
      status_ = LINEAR;
      setCachedLinearMatrix(cachedLinearFactor_in);
      // factorDim = cachedLinearFactor->dim();
      for(Key unmappedKey : cachedLinearFactor_in->keys()) {
        remappedKeys_.push_back(etree->getRemapKey(unmappedKey));
      }
    }
    else {
      // factorDim = nonlinearFactor_->dim();
      for(Key unmappedKey : nonlinearFactor_->keys()) {
        remappedKeys_.push_back(etree->getRemapKey(unmappedKey));
      }
    }

    // Add rhs to remapped keys
    remappedKeys_.push_back(0);

    updateBlockIndices();
    updateLowestHighestKeys();
  }

  void updateBlockIndices() {
    assert(status_ != REMOVING && status_ != REMOVED);
    blockIndices_.clear();
    assert(remappedKeys_.back() == 0);
    size_t curCol = 0;
    for(RemappedKey key : remappedKeys_) {
      size_t width = etree->colWidth(key);
      blockIndices_.push_back({key, curCol, width});
      curCol += width;
    }
  }

  size_t factorIndex() const { return factorIndex_; }

  const BlockIndexVector& blockIndices() const { return blockIndices_; }

  void updateLowestHighestKeys() {
    lowestKey_ = remappedKeys_.front();
    highestKey_ = remappedKeys_.front();
    for(auto it = remappedKeys_.begin(); it != remappedKeys_.end() - 1; it++) {
      // No need to check key 0
      RemappedKey key = *it;

      if(etree->orderingLess_(key, lowestKey_)) {
        lowestKey_ = key;
      }
      else if(etree->orderingLess_(highestKey_, key)) {
        highestKey_ = key;
      }
    }
  }

  RemappedKey lowestKey() { return lowestKey_; }
  RemappedKey highestKey() { return highestKey_; }

  const size_t dim() const { 
    /*return factorDim;*/ 
    throw std::runtime_error("dim() not implemented"); 
    return 0; 
  }

  // Copy the matrix in the linear factor into the cachedLinearMatrix
  // and set linearFactorType to be JACOBIAN or HESSIAN
  void setCachedLinearMatrix(sharedLinearFactor linearFactor);

  const FactorStatus& status() const { return status_; }

  void setStatusRemoved() {
    if(status_ == REMOVING || status_ == UNLINEARIZED) {
      status_ = REMOVED;
    }
    else {
      assert(0);
      throw std::runtime_error("Removing factor that is not REMOVING or UNLINEARIZED! Factor status");
    }
  }

  void setStatusRelinearize() {
    assert(status_ == LINEARIZED);
    status_ = RELINEARIZE;
  }

  void setStatusLinear() {
    assert(status_ == LINEARIZED || status_ == LINEAR || status_ == REMOVED);
    if(status_ != REMOVED) { status_ = LINEAR; }
    nonlinearFactor_ = nullptr;
  }

  void setStatusRemoving() {
    if(status_ != REMOVED) { status_ = REMOVING; }
  }

  bool hasMarginalizedKeys() const;

  // Linearize the factor if it needs to be linearized. 
  // Returns true if action is done, false otherwise
  // This does not change the status of the factor yet
  bool linearizeIfNeeded(const Values& theta) {

    assert(status_ != REMOVING && status_ != REMOVED);
    if(status_ == UNLINEARIZED || status_ == RELINEARIZE) {
      auto cachedLinearFactor = nonlinearFactor_->linearize(theta);
      setCachedLinearMatrix(cachedLinearFactor);

      status_ = LINEARIZED;
      
      return true;
    }
    return false;
  }

  sharedFactor nonlinearFactor() { return nonlinearFactor_; }

  LinearFactorType getLinearFactorType() {
    return linearFactorType_;
  }

  ColMajorMatrix<GEMMINI_TYPE>& getCachedMatrix() { return cachedLinearMatrix_; }

  void addNewKeys(const std::vector<RemappedKey>& newKeys) {
    assert(status_ != REMOVING && status_ != REMOVED);
    // We actually don't need to do a lot here as the nonlinearFactor will have already been 
    // updated with new keys, but we do need to relinearize the factor though
    status_ = RELINEARIZE;
    for(RemappedKey key : newKeys) {
      remappedKeys_.push_back(key);
    }

    // TODO: verify if new keys are added to the end of the list of keys
    throw std::runtime_error("factor add new keys");

    updateBlockIndices();
    updateLowestHighestKeys();

    return;
  }

  void markAffectedKeys(RemappedKeySet* affectedKeys);

  const RemappedKeyVector& remappedKeys() {
    return remappedKeys_;
  }

  void checkInvariant() const {}

  const KeyVector& keys();

  // Given the factor is linear or already linearized, ignore the blocks that correspond to 
  // keys that are marginalized out. Returns true if al keys are marginalized out
  // false otherwise
  bool marginalizeKeys();

  // const JacobianFactor* toJacobianFactor() const;

  // const HessianFactor* toHessianFactor() const;

  template<typename MATRIX, typename INFO, typename PREDICATE>
  void updateHessianJacobian(
      MATRIX& m, 
      double sign, 
      const INFO& info,
      const PREDICATE& pred) {
      
    // First check if any of the columns need to be processed
    bool skip = true;
    for(size_t i = 0; i < blockIndices_.size() - 1; i++) {
      // Higher key represents the column. Don't need last column
      const auto&[destR1, ordering1, status1] = info(i);
      // std::cout << "checking lowerKey = " << lowerKey << std::endl;
      if(pred(status1)) {
        skip = false;
        continue;
      }
    }
    if(skip) { return; }

    const auto& Ab = cachedLinearMatrix_;
    size_t height = Ab.rows();
    size_t width = Ab.cols();

    // Allocate a large scratch space
    std::vector<GEMMINI_TYPE> C(width * width, 0);

    // We want to get H^T = BB^T, where B = A^T
    syrk(width, width, height,
         &Ab(0, 0), &Ab(0, 0), C.data(),
         height, height, width,
         sign, 1, 
         false, true);

    size_t mrows = m.rows();
    size_t mcols = m.cols();
    
    for(size_t i = 0; i < blockIndices_.size() - 1; i++) {
      // Higher key represents the column. Don't need last column
      const auto&[lowerKey, srcCol1, srcW1] = blockIndices_.at(i);
      const auto&[destR1, ordering1, status1] = info(i);
      if(!pred(status1)) {
        continue;
      }
      for(size_t j = 0; j < blockIndices_.size(); j++) {
        const auto&[higherKey, srcCol2, srcW2] = blockIndices_.at(j);
        const auto&[destR2, ordering2, status2] = info(j);
        if(ordering2 < ordering1) {
          continue;
        }

        Eigen::Block<MATRIX> destBlock(m, destR2, destR1, srcW2, srcW1);

        if(srcCol2 >= srcCol1) {
          if(destR2 >= destR1) {
            add_block(C.data(), &m(0, 0), 
                      width, mrows, 
                      srcCol2, srcCol1, 
                      destR2, destR1, 
                      srcW2, srcW1);
          }
          else {
            add_block_transpose(C.data(), &m(0, 0), 
                                width, mrows, 
                                srcCol2, srcCol1, 
                                destR1, destR2, 
                                srcW2, srcW1);
          }
        }
        else {
          if(destR2 >= destR1) {
            add_block_transpose(C.data(), &m(0, 0), 
                                width, mrows, 
                                srcCol1, srcCol2, 
                                destR2, destR1, 
                                srcW1, srcW2);
          }
          else {
            add_block(C.data(), &m(0, 0), 
                      width, mrows, 
                      srcCol1, srcCol2, 
                      destR1, destR2, 
                      srcW1, srcW2);
          }
        }
      }
    }

  }

  template<typename MATRIX, typename INFO, typename PREDICATE>
  void updateHessianHessian(
      MATRIX& m, 
      double sign, 
      const INFO& info,
      const PREDICATE& pred) {
    const auto& AbtAb = cachedLinearMatrix_;
    for(size_t i = 0; i < blockIndices_.size() - 1; i++) {
      // Higher key represents the column. Don't need last column
      const auto&[lowerKey, srcCol1, srcW1] = blockIndices_.at(i);
      const auto&[destR1, ordering1, status1] = info(i);
      // std::cout << "checking lowerKey = " << lowerKey << std::endl;
      if(!pred(status1)) {
        continue;
      }
      for(size_t j = 0; j < blockIndices_.size(); j++) {
        const auto&[higherKey, srcCol2, srcW2] = blockIndices_.at(j);
        const auto&[destR2, ordering2, status2] = info(j);
        if(ordering2 < ordering1) {
          assert(destR2 < destR1);
          continue;
        }

        Eigen::Block<MATRIX> destBlock(m, destR2, destR1, srcW2, srcW1);

        assert(destR2 >= destR1);

        if(srcCol2 >= srcCol1) {
          // Note: We need to access the upper triangular part of the Hessian matrix
          Eigen::Block<const GemminiMatrix> AbtAb_ij(AbtAb, srcCol1, srcCol2, srcW1, srcW2);

          destBlock.noalias() += sign * AbtAb_ij.transpose();
        }
        else {
          // Note: We need to access the upper triangular part of the Hessian matrix
          Eigen::Block<const GemminiMatrix> AbtAb_ij(AbtAb, srcCol2, srcCol1, srcW2, srcW1);

          destBlock.noalias() += sign * AbtAb_ij;
        }

      }
    }

  }

  // Populate the columns of the MATRIX m with the sign * Hessian of the factor
  // skipping column key if pred(key) == false
  template<typename MATRIX, typename INFO, typename PREDICATE>
  void updateHessian(
      MATRIX& m, 
      double sign, 
      const INFO& info,
      const PREDICATE& pred) {

    // std::cout << "Updating hessian for " << *this << " sign = " << sign << std::endl;

    if(linearFactorType_ == JACOBIAN) {
      updateHessianJacobian(m, sign, info, pred);
    }
    else if(linearFactorType_ == HESSIAN) {
      updateHessianHessian(m, sign, info, pred);   
    }
    else {
      throw std::runtime_error("Factor cannot be dynamically cast");
    }

  }


  void printKeys(std::ostream& os) const;
};

} // namespace gtsam
