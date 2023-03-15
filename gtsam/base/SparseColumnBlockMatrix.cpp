/**
* @file    SparseColumnBlockMatrix.cpp
* @brief   Represents a column matrix of fixed width and predefined block heights
* @author  Roger Hsiao (rogerhh)
* @date    Feb. 8, 2023
*/

#include <cassert>
#include <gtsam/base/SparseColumnBlockMatrix.h>

using namespace std;

namespace gtsam {

using Block = SparseColumnBlockMatrix::Block;
using constBlock = SparseColumnBlockMatrix::constBlock;
using RowHeightPair = SparseColumnBlockMatrix::RowHeightPair;

const size_t SparseColumnBlockMatrix::width() const {
    return width_; 
}

const size_t SparseColumnBlockMatrix::height() const {
    assert(newMaxHeight_ == maxHeight_);
    return newMaxHeight_;
}

const Key SparseColumnBlockMatrix::key() const {
   return key_;
}

size_t SparseColumnBlockMatrix::matrixHeight() const {
    return maxHeight_;
}

// if underlying matrix height is the same as expected matrix height
bool SparseColumnBlockMatrix::allocated() const {
    return maxHeight_ == newMaxHeight_;
}

// Diagonal flag to indicate if diagonal block should be initialized
// All column matrices should have a diagonal block corresponding to the column index
// In the Cholesky matrix, we want the diagonal block to be on top
// In the Hessian matrix, we don't care about the ordering of blocks
// So we should always keep diagonal block on top
SparseColumnBlockMatrix::SparseColumnBlockMatrix(
        const Key key_in, const size_t width_in, bool diagonal) 
        : key_(key_in), width_(width_in) {
    if(diagonal) {
        preallocateBlock(key_, width_, true);
    }
}

// bool SparseColumnBlockMatrix::tryAllocateBlock(const Key otherKey,
//                                                const size_t height) {
//     auto iterPair = blockStartMap_.insert({otherKey, {newMaxHeight_, height}});
//     if(iterPair.second) {
//         // If inserted, increase max height and set vector 
//         newMaxHeight_ += height;
//         blockStartVec_.push_back({otherKey, height});
//     }
//     return iterPair.second;
// }

bool SparseColumnBlockMatrix::preallocateBlock(const Key otherKey,
                                               const size_t height,
                                               const bool initialize) {
    auto iterPair = blockStartMap_.insert({otherKey, {newMaxHeight_, height}});
    if(iterPair.second) {
        // If inserted, increase max height and set vector 
        blockStartVec_.push_back({otherKey, {newMaxHeight_, height}});
        newMaxHeight_ += height;
    }
    else if(initialize) {
        // If did not insert and initialize, set block to 0
        if(iterPair.first->second.first + iterPair.first->second.second <= maxHeight_) {
            blockRange(iterPair.first->second.first, iterPair.first->second.second).setZero();
        }
    }
    return iterPair.second;
}

// Actually allocate and initialize the amount we need
void SparseColumnBlockMatrix::resolveAllocate() {
    assert(newMaxHeight_ >= maxHeight_);    
    if(newMaxHeight_ == maxHeight_) {
        return;
    }
    if(maxHeight_ == 0) {
        matrix_ = RowMajorMatrix(newMaxHeight_, width_);
    }
    else {
        matrix_.conservativeResize(newMaxHeight_, Eigen::NoChange_t());
    }
    matrix_.block(maxHeight_, 0, newMaxHeight_ - maxHeight_, width_).setZero();
    maxHeight_ = newMaxHeight_;
}

void SparseColumnBlockMatrix::setZero() {
    matrix_.setZero();
}

void SparseColumnBlockMatrix::setZero(const Key i) {
    block(i).setZero();
}

// reset blockStart* assignments except for the diagonal block, 
// but don't touch the underlying matrix in case we need that
// memory later. Excess memory will be freed by resize
void SparseColumnBlockMatrix::resetBlocks(bool diagonal) {
    maxHeight_ = 0;
    newMaxHeight_ = 0;
    blockStartVec_.clear();
    blockStartMap_.clear();
    if(diagonal) {
        preallocateBlock(key_, width_, true);
    }
}

bool SparseColumnBlockMatrix::blockExists(const Key i) const {
    return blockStartMap_.find(i) != blockStartMap_.end();
}

// Access functions
Block SparseColumnBlockMatrix::block(const Key i) {
    auto pair = blockStartMap_.at(i); 
    return matrix_.block(pair.first, 0, pair.second, width_);
}

const constBlock SparseColumnBlockMatrix::block(const Key i) const {
    auto pair = blockStartMap_.at(i); 
    return matrix_.block(pair.first, 0, pair.second, width_);
}

// Access the underlying matrix with variable height but fixed width
// This is used for when we want to compute the contribution blocks
Block SparseColumnBlockMatrix::blockRange(
        const size_t startRow, const size_t height) {
    return matrix_.block(startRow, 0, height, width_);
}

const constBlock SparseColumnBlockMatrix::blockRange(
        const size_t startRow, const size_t height) const {
    return matrix_.block(startRow, 0, height, width_);
}

Block SparseColumnBlockMatrix::submatrix(const size_t startRow, const size_t startCol,
                                         const size_t height, const size_t width) {
    return matrix_.block(startRow, startCol, height, width);
}

Block SparseColumnBlockMatrix::diagonalBlock() {
    assert(maxHeight_ >= width_);
    return matrix_.block(0, 0, width_, width_);
}
const constBlock SparseColumnBlockMatrix::diagonalBlock() const {
    assert(maxHeight_ >= width_);
    return matrix_.block(0, 0, width_, width_);
}

bool SparseColumnBlockMatrix::hasBelowDiagonalBlocks() const {
    return newMaxHeight_ > width_; // if matrix doesn't just have the diagonal blocks
}

Block SparseColumnBlockMatrix::belowDiagonalBlocks() {
    return matrix_.block(width_, 0, maxHeight_ - width_, width_);
}
const constBlock SparseColumnBlockMatrix::belowDiagonalBlocks() const {
    return matrix_.block(width_, 0, maxHeight_ - width_, width_);
}

const vector<pair<Key, RowHeightPair>>& SparseColumnBlockMatrix::blockStartVec() const {
    return blockStartVec_;
}
const unordered_map<Key, RowHeightPair>& SparseColumnBlockMatrix::blockStartMap() const {
    return blockStartMap_;
}

void SparseColumnBlockMatrix::print(ostream& os) const {
    for(const auto p : blockStartVec_) {
        Key k = p.first;
        size_t row = p.second.first;
        size_t height = p.second.second;
        os << "Key: " << k << "\n";
        os << blockRange(row, height) << "\n\n";
    }
}

void SparseColumnBlockMatrix::checkInvariant() const {
    assert(blockStartVec_.size() == blockStartMap_.size());
    for(auto p : blockStartVec_) {
        Key k = p.first;
        assert(blockStartMap_.find(k) != blockStartMap_.end());
        assert(blockStartMap_.at(k) == p.second);
    }
}

} // namespace gtsam
