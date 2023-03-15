/**
* @file    SparseVerticalBlockMatrix.h
* @brief   Access to matrices via blocks of pre-defined sized. Used as a unified Jacobian matrix
* @author  Roger Hsiao (rogerhh)
* @date    Feb. 8, 2023
*/
#pragma once

#include <gtsam/base/Matrix.h>
#include <gtsam/base/MatrixSerialization.h>
#include <gtsam/base/FastVector.h>
#include <Eigen/Core>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "base/BlockMatrix.h"

namespace gtsam {

  /**
   * This class stores a sparse block matrix and allows it to be accessed as a collection of vertical
   * blocks. More blocks can be added to each column later
   *
   * @ingroup base */
    class GTSAM_EXPORT JacobianMatrix
    {
    public:
        // Allocate the Jacobian matrix for a factor
        // takes in a vector of indices in the factor
        // A factor index (row number J)
        // number of rows in that factor
      void allocateJacobianFactor(const std::vector<Key>& varIndices, 
                                  const size_t factorIndex, 
                                  const size_t numRows);

      class RowIterator {};

      class ColIterator {};

      RowIterator begin();
      RowIterator end();
      const RowIterator begin();
      const RowIterator end();
      ColIterator begin();
      ColIterator end();
      const ColIterator begin();
      const ColIterator end();
      



    private:
      // This vector is column major. values_[key] gives the column vector of some variable
      std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> columnMatrices_; 
      // Each column is represented by a Matrix object. Use RowMajor Storage to support growing the matrix
      std::vector<size_t> colWidths_; // Use column width = 0 to indicate unused/deleted variables
      std::vector<std::vector<size_t>> rowIndices_; // Each column needs to store the row indices
      std::vector<std::vector<size_t>> factorIndices_; // Factor indices of each column
      std::vector<std::unordered_map<size_t, size_t>> factorToRowIndices_;
      // Use Matrix.conservativeResize(NoChange_t, rows) to resize matrices




      // // Allocate num_rows * column_widths_[i] doubles for each index
      // void allocateJacobianFactor(const std::vector<Key>& var_indices, 
      //                             const size_t factor_index, 
      //                             const size_t num_rows) {
      //     for(const Key key : var_indices) {
      //         blocks_[key][factor] = values_[key].size();
      //         values_[key].resize(values_[key].size() + num_rows * col_widths_[key]);
      //     }
      //     row_widths_.push_back(num_rows);
      // }

      // Returns view of block at location {Key, factor_index} 
      BlockMatrix block(const Key key, const size_t factor_index) {
        return blocks_[key][factor_index];
      }


  };
}
