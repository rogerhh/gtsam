/**
* @file    SparseAdjointBlockMatrix.h
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

namespace gtsam {

  /**
   * This class stores a sparse adjoint block matrix and allows it to be accessed as a collection of blocks columns
   * Can be used as a symmetric matrix
   * Blocks are accesses by two keys i, j, and only block[i][j] where i <= j will be populated. 
   * The caller needs ensure that i <= j and transpose the returned block if necessary
   * No selfadjointView on diagonal needed because the matrices we're dealing with are so small
   *
   * @ingroup base */
    
class GTSAM_EXPORT SparseAdjointBlockMatrix {
private:
    typedef SparseAdjointBlockMatrix This;

    // The underlying storage of an adjoint matrix is a collection of column blocks
    std::vector<SparseColumnBlockMatrix> columns_;

public:
    using SparseColumnBlockMatrix::Block;
    using SparseColumnBlockMatrix::constBlock;

    void addColumn(const Key key, const size_t width) {
    }

    void preallocateBlock(const size_t i, const size_t j) {
        if(i <= j) {
            upperTriangular.preallocateOrInitialize(i, j, initialize);
        }
        else {
            upperTriangular.preallocateOrInitialize(j, i, initialize);
        }
    }

    void resolveAllocate() {
    }

    void resetColumn(const Key key) {
        upperTriangular.resetColumn(key);
    }

    Block block(const Key i, const Key j) {
        if(i <= j) {
            return upperTriangular.block(i, j);
        }
        return upperTriangular.block(j, i);
    }

    const constBlock block(const Key i, const Key j) const {
        if(i <= j) {
            return upperTriangular.block(i, j);
        }
        return upperTriangular.block(j, i);
    }

};

}
