/**
* @file    LowerTriangularColumnMatrix.h
* @brief   Represents a column for a lower triangular matrix, contains a diagonal block and a row at the end for Atb. Column blocks need to support insertion in the middle. We don't support it for now
* @author  Roger Hsiao (rogerhh)
* @date    Mar. 15, 2023
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

namespace gtsam {

class GTSAM_EXPORT LowerTriangularColumnMatrix {
public:
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorMatrix;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> ColMajorMatrix;
    typedef Eigen::Block<ColMajorMatrix> Block;
    typedef Eigen::Block<const ColMajorMatrix> constBlock;
    
    // Each key maps to {blockStartRow, blockHeight}
    typedef std::pair<size_t, size_t> RowHeightPair;

    typedef size_t RowKey; // -1 indicates the Atb row

    typedef std::vector<std::pair<RowKey, RowHeightPair>> BlockIndexVector;

    static const size_t LAST_ROW;

private:

    Key key_;
    size_t width_ = 0;   // Use 0 to indicate invalid
    size_t maxHeight_ = 0;
    size_t newMaxHeight_ = 0;   // Height not including last row

    // newMaxHeight_ is always height of the matrix - 1
    // vector for fast sequential access. This is mostly used in the Cholesky block
    // where the ordering is fixed and there should be no empty blocks
    BlockIndexVector blockIndices_;

    // Underlying storage. Row major to support growing rows
    ColMajorMatrix matrix_;

public:
    
    // Atb row will always be populated
    LowerTriangularColumnMatrix(const RowKey key_in, const size_t width_in);

    // Access parameter functions
    const size_t width() const;

    // This returns the nominal height of the matrix
    const size_t height() const;
    
    // Height of the underlying matrix. If 0, matrix is unallocated
    size_t matrixHeight() const;

    const Key key() const;
    
    // if underlying matrix height is the same as expected matrix height
    bool allocated() const;

    void preallocateBlocks(const BlockIndexVector& blocks, size_t startIndex = 0);
    
    // Add new blocks to the end of the column
    void preallocateBlocks(const std::vector<std::pair<RowKey, size_t>>& blocks, size_t startIndex = 0);

    // Add new block to the end of the column, new block must not already be in column
    void preallocateBlock(const std::pair<RowKey, size_t>& blocks);

    // Actually allocate and initialize the amount we need
    void resolveAllocate();

    // Don't change block structure, just set everything to 0
    void setZero();

    // reset blockStart* assignments except for the diagonal block, 
    // but don't touch the underlying matrix in case we need that
    // memory later. Excess memory will be freed by resize
    void resetBlocks();

    // Access the underlying matrix with variable height
    // This is used for when we want to compute the contribution blocks
    Block lastRow();
    Block blockRange(const size_t startRow, const size_t height);
    const constBlock blockRange(const size_t startRow, const size_t height) const;

    const BlockIndexVector& blockIndices() const; 

    void print(std::ostream& os) const;
    void printBlockIndices(std::ostream& os) const;

    void checkInvariant() const;

    // Reorder blocks according to reorderedKeys. reorderedKeys should
    // only contain keys that need to be reordered
    // oldLowestKeyIndex is the position of the lowest key
    void reorderBlocks(const std::vector<Key>& reorderedKeys, size_t oldLowestKeyIndex);

private:

    void preallocateLastRow();
}; 


} // namespace gtsam
