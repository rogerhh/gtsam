/**
* @file    SparseSymmetricBlockMatrix.h
* @brief   Access to matrices via blocks of pre-defined sized. Used as a unified Hessian matrix
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
#include <iostream>

namespace gtsam {

  /**
   * This class stores a sparse symmetric block matrix and allows it to be accessed as a collection of blocks columns
   * Blocks are accesses by two keys i, j, and only block[i][j] where i >= j will be populated. 
   * Effectively storing the lower triangular block
   * The caller needs make sure to transpose the returned block based on
   * No selfadjointView needed because the matrices we're dealing with are so small
   *
   * @ingroup base */

class GTSAM_EXPORT SparseSymmetricBlockMatrix {
public:
    typedef SparseColumnBlockMatrix::Block Block;
    typedef SparseColumnBlockMatrix::constBlock constBlock;
    typedef SparseColumnBlockMatrix::RowHeightPair RowHeightPair;

private:
    typedef SparseSymmetricBlockMatrix This;
    std::vector<SparseColumnBlockMatrix> columns_;
    std::unordered_set<size_t> needAlloc_;

public:
    void addColumn(const Key key, const size_t width);

    size_t colWidth(const Key key) const;
    size_t colHeight(const Key key) const;

    void preallocateBlock(const Key i, const Key j, bool initialize=true);

    void resolveAllocate();

    void setZero(const Key i, const Key j);
    void setZero(const Key col);

    void resetColumn(const Key col);

    bool blockExists(const Key i, const Key j) const;

    Block block(const Key i, const Key j);
    const constBlock block(const Key i, const Key j) const;

    Block colBlockRange(const Key key, const size_t row, const size_t height);
    const constBlock colBlockRange(
            const Key key, const size_t row, const size_t height) const;

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

}
