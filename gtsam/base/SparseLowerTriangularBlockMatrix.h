/**
* @file    SparseLowerTriangularBlockMatrix.h
* @brief   Access to matrices via blocks of pre-defined sized. Used as a unified Cholesky factor matrix or a underlying storage to a sparse symmetric block matrix
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
#include <gtsam/base/SparseColumnBlockMatrix.h>

namespace gtsam {

  /**
   * This class stores a sparse lower triangular block matrix and allows it to be accessed as a collection of column blocks.
   * There is no check on the indices because the columns can be reordered. However, the user must make sure that (i, j) and (j, i) are not nonzero at the same time
   *
   * @ingroup base */

class GTSAM_EXPORT SparseLowerTriangularBlockMatrix {
public:
    typedef SparseColumnBlockMatrix::Block Block;
    typedef SparseColumnBlockMatrix::constBlock constBlock;
    typedef SparseColumnBlockMatrix::RowHeightPair RowHeightPair;

private:
    typedef SparseLowerTriangularBlockMatrix This;
    std::vector<SparseColumnBlockMatrix> columns_;
    std::unordered_set<size_t> needAlloc_;

public:
    void addColumn(const Key key, const size_t width);

    size_t colWidth(const Key key) const;
    size_t colHeight(const Key key) const;

    bool preallocateBlock(const Key i, const Key j, bool initialize=true);

    void resolveAllocate();
    void resolveAllocate(const Key key);

    void setZero(const Key i, const Key j);
    void setZero(const Key col);

    void resetColumn(const Key col);

    bool blockExists(const Key i, const Key j) const;

    Block block(const Key i, const Key j);
    const constBlock block(const Key i, const Key j) const;

    Block colBlockRange(const Key key, const size_t row, const size_t height);
    const constBlock colBlockRange(
            const Key key, const size_t row, const size_t height) const;

    Block colDiagonalBlock(const Key key);
    const constBlock colDiagonalBlock(const Key key) const;

    Block colBelowDiagonalBlocks(const Key key);
    const constBlock colBelowDiagonalBlocks(const Key key) const;

    SparseColumnBlockMatrix& column(const Key key);

    void print(std::ostream& os) const;

    void printColumn(std::ostream& os, const Key i) const;

    void checkInvariant() const;

    // void addColumn(const Key key, const size_t width) {
    //     upperTriangular.addColumn(key, width);
    // }

    // void preallocateOrInitialize(const Key i, const Key j, const bool initialize) {
    //     if(i <= j) {
    //         upperTriangular.preallocateOrInitialize(i, j, initialize);
    //     }
    //     else {
    //         upperTriangular.preallocateOrInitialize(j, i, initialize);
    //     }
    // }

    // void resolveAllocate() {
    //     upperTriangular.resolveAllocate();
    // }

    // void resetColumn(const Key key) {
    //     upperTriangular.resetColumn(key);
    // }

    // Block block(const Key i, const Key j) {
    //     if(i <= j) {
    //         return upperTriangular.block(i, j);
    //     }
    //     return upperTriangular.block(j, i);
    // }

    // const constBlock block(const Key i, const Key j) const {
    //     if(i <= j) {
    //         return upperTriangular.block(i, j);
    //     }
    //     return upperTriangular.block(j, i);
    // }

};

// class GTSAM_EXPORT SparseLowerTriangularBlockMatrix {
// public:
//     typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorMatrix;
//     typedef Eigen::Block<RowMajorMatrix> Block;
//     typedef Eigen::Block<const RowMajorMatrix> constBlock;
// 
// private:
// 
//     struct ColumnMatrix {
// 
//         Key key;
//         size_t width = 0;   // Use 0 to indicate invalid
//         size_t maxHeight = 0;
//         size_t newMaxHeight = 0;
//         std::unordered_map<Key, std::pair<size_t, size_t>> blockStart;  // Each key maps to {blockStartRow, height}
//         RowMajorMatrix matrix;  // Row major to support growing
// 
//         // TODO: Initialize a column matrix with a diagonal block
//         ColumnMatrix(const Key key_in, const size_t width_in);
// 
//         // Check if block exists. If not, set appropriate indices to be allocated later
//         // returns true if allocated
//         bool preallocateOrInitialize(const Key otherKey, 
//                                      const size_t height,
//                                      const bool initialize);
// 
//         // Actually allocate and initialize the amount we need
//         // In case of relinearization, maxHeight = 0. We might free up memory by resizing
//         void resolveAllocate();
// 
//         // reset block start assignments, but don't touch the underlying matrix
//         // in case we need that memory later. Excess memory will be freed by
//         // conservativeResize
//         void resetBlocks();
// 
//         Block block(const Key i);
// 
//         const constBlock block(const Key i) const;
//     };
// 
// public:
//     void addColumn(const Key key, const size_t width);
// 
//     void preallocateOrInitialize(const Key i, const Key j, const bool initialize);
// 
//     void resolveAllocate();
// 
//     void resetColumn(const Key key);
// 
//     Block block(const Key i, const Key j);
// 
//     const constBlock block(const Key i, const Key j) const;
// 
// private:
//     // Each column is represented by a Matrix object. Use RowMajor Storage to support growing the matrix
//     std::vector<ColumnMatrix> columnMatrices_; 
//     std::unordered_set<size_t> needAlloc;
// 
// };

} // namespace gtsam
