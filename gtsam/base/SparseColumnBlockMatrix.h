/**
* @file    SparseColumnBlockMatrix
* @brief   Represents a column matrix of fixed width and predefined block heights. No predefined diagonal block
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

namespace gtsam {

  /**
   * This class stores a vertical block matrix
   *
   * @ingroup base */

class GTSAM_EXPORT SparseColumnBlockMatrix {
public:
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMajorMatrix;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> ColMajorMatrix;
    typedef Eigen::Block<ColMajorMatrix> Block;
    typedef Eigen::Block<const ColMajorMatrix> constBlock;
    
    // Each key maps to {blockStartRow, blockHeight}
    typedef std::pair<size_t, size_t> RowHeightPair;

private:

    typedef size_t RowKey;

    Key key_;
    size_t width_ = 0;   // Use 0 to indicate invalid
    size_t maxHeight_ = 0;
    size_t newMaxHeight_ = 0;
    // vector for fast sequential access. This is mostly used in the Cholesky block
    // where the ordering is fixed and there should be no empty blocks
    std::vector<std::pair<RowKey, RowHeightPair>> blockStartVec_;

    // Unordered map for random access. This is mostly used in the Hessian block
    // where is ordering is fluid and there may be empty blocks 
    std::unordered_map<RowKey, RowHeightPair> blockStartMap_;

    // Underlying storage. Row major to support growing rows
    ColMajorMatrix matrix_;

public:

    // Access parameter functions
    const size_t width() const;

    // This returns the nominal height of the matrix
    const size_t height() const;

    const Key key() const;
    
    // In the Hessian matrix, we don't care about the ordering of blocks
    SparseColumnBlockMatrix(const Key key_in, const size_t width_in);

    // Check if blocks exist. If not, set appropriate indices to be allocated/initialized later
    // returns true of allocated
    // Requires blocks to be sorted
    // We don't want to allocate now because we might need to move data around
    bool preallocateBlock(const Key otherKey,
                          const size_t height,
                          const bool initialize);

    // Actually allocate and initialize the amount we need
    void resolveAllocate();

    // Don't change block structure, just set everything to 0
    void setZero();

    // Set a single block to 0
    void setZero(const Key i);

    // reset blockStart* assignments 
    // but don't touch the underlying matrix in case we need that
    // memory later. Excess memory will be freed by resize
    void resetBlocks();

    bool blockExists(const Key i) const;

    // Access functions
    Block block(const Key i);
    const constBlock block(const Key i) const;

    // Access the underlying matrix with variable height
    // This is used for when we want to compute the contribution blocks
    Block blockRange(const size_t startRow, const size_t height);
    const constBlock blockRange(const size_t startRow, const size_t height) const;

    const std::vector<std::pair<Key, RowHeightPair>>& blockStartVec() const; 
    const std::unordered_map<Key, RowHeightPair>& blockStartMap() const; 

    void print(std::ostream& os) const;

    void checkInvariant() const;

    // Height of the underlying matrix. If 0, matrix is unallocated
    size_t matrixHeight() const;
    
    // if underlying matrix height is the same as expected matrix height
    bool allocated() const;

    // Reorder blocks according to reorderedKeys. reorderedKeys should
    // only contain keys that need to be reordered
    // oldLowestKeyIndex is the position of the lowest key
    void reorderBlocks(const std::vector<RowKey>& reorderedKeys, 
                       size_t oldLowestKeyIndex);


private:
    // Returns true if block did not exist
    // bool tryAllocateBlock(const Key otherKey, const size_t height);

}; 

} // namespace gtsam
