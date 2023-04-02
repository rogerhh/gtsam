/**
* @file    JacobianMatrix
* @brief   Represents a collection of column matrices
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
#include <iostream>
#include <gtsam/base/SparseColumnBlockMatrix.h>

namespace gtsam {

  /**
   * This class stores collection of vertical block matrices
   *
   * @ingroup base */

class GTSAM_EXPORT JacobianMatrix {
private:
    std::vector<SparseColumnBlockMatrix> columns_;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorMatrix;

public:
    typedef SparseColumnBlockMatrix::Block Block;
    typedef SparseColumnBlockMatrix::constBlock constBlock;

    JacobianMatrix();

    void addColumn(const Key key, const size_t width);

    // FactorIndex is the row, and Key is the column
    void preallocateBlock(const FactorIndex i, const Key j, const size_t height, bool initialize);

    void resolveAllocate(const Key i);

    void setZero(const Key i);

    void setZero(const Key i, const Key j);

    void resetColumn(const Key i);

    bool blockExists(const Key i, const Key j) const;

    // Access functions
    Block block(const Key i, const Key j);
    const constBlock block(const Key i, const Key j) const;

    // Access the underlying matrix with variable height
    // This is used for when we want to compute the contribution blocks
    Block blockRange(const size_t startRow, const size_t height);
    const constBlock blockRange(const size_t startRow, const size_t height) const;

    void print(std::ostream& os);

};

} // namespace gtsam
