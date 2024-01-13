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

namespace gtsam {
    
  typedef Key RemappedKey;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> ColMajorMatrix;
  typedef Eigen::VectorXd Vector;
  typedef std::tuple<RemappedKey, size_t, size_t> KeyRowHeight;
  typedef std::tuple<size_t, size_t> RowHeight;
  typedef std::vector<KeyRowHeight> BlockIndexVector;
  typedef std::unordered_map<RemappedKey, RowHeight> BlockIndexMap;
  typedef std::vector<RemappedKey> RemappedKeyVector;
  typedef std::set<RemappedKey> RemappedKeySet;

  enum MarkedStatus {UNMARKED=0, EDIT=1, RECONSTRUCT=2, NEW=3, MARGINALIZED=4};
  enum FactorStatus {UNLINEARIZED=0, LINEARIZED=1, RELINEARIZE=2, LINEAR=3, 
                     REMOVING=4, REMOVED=5, UNUSED=6};

  static constexpr size_t BLOCK_INDEX_KEY = 0;
  static constexpr size_t BLOCK_INDEX_ROW = 1;
  static constexpr size_t BLOCK_INDEX_HEIGHT = 2;

} // namespace gtsam
