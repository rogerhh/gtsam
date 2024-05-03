/**
* @file    Workspace.h
* @brief   Work space memory for multifrontal elimination
* @author  Roger Hsiao (rogerhh)
* @date    Apr. 12, 2023
*/

#pragma once

#include <Eigen/Dense>
#include <Eigen/Core>
#include <cassert>
#include <vector>
#include <iostream>
#include <utility>
#include <tuple>
#include <gtsam/linear/CholeskyEliminationTreeTypes.h>

#include <gtsam/linear/gemmini_functions.h>

namespace gtsam {

class Workspace {
public:
  // typedef std::tuple<size_t, size_t, size_t> BlockIndexVector;
  // typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> ColMajorMatrix;

private:

  GEMMINI_TYPE* data = nullptr;
  GEMMINI_TYPE* stack_ptr = nullptr;
  size_t size = 0;

  // std::vector<double*> stack_ptrs;
  std::vector<std::vector<GEMMINI_TYPE*>> stack_ptrs;
  std::vector<std::vector<Eigen::Map<ColMajorMatrix<GEMMINI_TYPE>>>> stack;

public:

  Workspace() {}

  ~Workspace() {
    deallocate();
  }

  void deallocate() {
    assert(stack_ptrs.empty());
    assert(stack.empty());
    free(data);
    data = nullptr;
    stack_ptr = nullptr;
    size = 0;
  }

  // allocate memory for underlying memory
  // This should only be run once after symbolic elim
  void allocate(size_t size_in) {
    assert(stack_ptrs.empty());
    data = (GEMMINI_TYPE*) realloc(data, size_in * sizeof(GEMMINI_TYPE));
    stack_ptr = data;
    size = size_in;
  }

  size_t push_matrices(const std::vector<std::pair<size_t, size_t>>& matrixSizes) {
    stack.push_back(std::vector<Eigen::Map<ColMajorMatrix<GEMMINI_TYPE>>>());
    stack_ptrs.push_back(std::vector<GEMMINI_TYPE*>());
    // stack_ptrs.push_back(stack_ptr);
    for(const auto& p : matrixSizes) {
      stack_ptrs.back().push_back(stack_ptr);
      size_t r, c;
      std::tie(r, c) = p;
      assert(stack_ptr + r * c <= data + size);
      stack.back().push_back(Eigen::Map<ColMajorMatrix<GEMMINI_TYPE>>(stack_ptr, r, c));
      stack.back().back().setZero();
      stack_ptr += r * c;
    }
    return stack.size() - 1;
  }

  const std::vector<Eigen::Map<ColMajorMatrix<GEMMINI_TYPE>>>& get_matrices(size_t index) {
    return stack.at(index);
  }

  const std::vector<GEMMINI_TYPE*>& get_ptrs(size_t index) {
    return stack_ptrs.at(index);
  }

  void pop_matrices() {
    stack_ptr = stack_ptrs.back().front();
    stack.pop_back();
    stack_ptrs.pop_back();
  }

  bool empty() const { return stack.empty(); }
};

}   // namespace gtsam
