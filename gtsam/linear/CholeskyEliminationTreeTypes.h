/**
* @file    CholeskyEliminationTreeTypes.h
* @brief   Common types
* @author  Roger Hsiao (rogerhh)
* @date    Jul. 10, 2023
*/

#pragma once

#include <set>
#include <unordered_map>
#include <gtsam/base/Matrix.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>

#include <gtsam/linear/gemmini_functions.h>

namespace gtsam {
    
  typedef Key RemappedKey;
  // typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> ColMajorMatrix;
  template <typename T=double>
  using ColMajorMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
  typedef Eigen::VectorXd Vector;
  typedef std::tuple<RemappedKey, size_t, size_t> KeyRowHeight;
  typedef std::tuple<size_t, size_t> RowHeight;
  typedef std::vector<KeyRowHeight> BlockIndexVector;
  typedef std::unordered_map<RemappedKey, RowHeight> BlockIndexMap;
  typedef std::vector<RemappedKey> RemappedKeyVector;
  typedef std::set<RemappedKey> RemappedKeySet;

  enum MarkedStatus {UNMARKED=0, EDIT=1, RECONSTRUCT=2, NEW=3, MARGINALIZED=4, FIXED=5};
  enum FactorStatus {UNLINEARIZED=0, LINEARIZED=1, RELINEARIZE=2, LINEAR=3, 
                     REMOVING=4, REMOVED=5, UNUSED=6};
  enum CostStatus {COST_UNMARKED = 0, COST_FIXED = 1, COST_MARKED = 2, COST_RELIN = 3};

  static constexpr size_t BLOCK_INDEX_KEY = 0;
  static constexpr size_t BLOCK_INDEX_ROW = 1;
  static constexpr size_t BLOCK_INDEX_HEIGHT = 2;

  typedef Eigen::Matrix<GEMMINI_TYPE, Eigen::Dynamic, Eigen::Dynamic> GemminiMatrix;
  typedef Eigen::Matrix<GEMMINI_TYPE, Eigen::Dynamic, 1> GemminiVector;

} // namespace gtsam
