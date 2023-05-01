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
#include <gtsam/base/LowerTriangularColumnMatrix.h>

namespace gtsam {

  /**
   * This class stores a sparse lower triangular block matrix and allows it to be accessed as a collection of column blocks.
   * There is no check on the indices because the columns can be reordered. However, the user must make sure that (i, j) and (j, i) are not nonzero at the same time
   *
   * @ingroup base */

class GTSAM_EXPORT SparseLowerTriangularBlockMatrix {
public:
    typedef LowerTriangularColumnMatrix::Block Block;
    typedef LowerTriangularColumnMatrix::constBlock constBlock;
    typedef LowerTriangularColumnMatrix::RowHeightPair RowHeightPair;
    static const size_t LAST_ROW;

private:
    typedef SparseLowerTriangularBlockMatrix This;
    std::vector<LowerTriangularColumnMatrix> columns_;

public:
    void addColumn(const Key key, const size_t width);

    size_t colWidth(const Key key) const;
    size_t colHeight(const Key key) const;

    void resolveAllocate();
    void resolveAllocate(const Key key);

    void setZero(const Key col);

    void resetColumn(const Key col);

    Block colBlockRange(const Key key, const size_t row, const size_t height);
    const constBlock colBlockRange(
            const Key key, const size_t row, const size_t height) const;

    LowerTriangularColumnMatrix& column(const Key key);

    LowerTriangularColumnMatrix& at(const Key key);
    const LowerTriangularColumnMatrix& at(const Key key) const;

    void print(std::ostream& os) const;
    void printColumn(std::ostream& os, const Key i) const;

    void checkInvariant() const;
};

} // namespace gtsam
