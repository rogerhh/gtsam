/**
* @file    LowerTriangularColumnMatrix.cpp
* @brief   Represents a column for a lower triangular matrix, contains a diagonal block and a row at the end for Atb. Column blocks need to support insertion in the middle
* @author  Roger Hsiao (rogerhh)
* @date    Mar. 15, 2023
*/

#include <cassert>
#include <gtsam/base/LowerTriangularColumnMatrix.h>

using namespace std;

namespace gtsam {

using Block = LowerTriangularColumnMatrix::Block;
using constBlock = LowerTriangularColumnMatrix::constBlock;
using RowHeightPair = LowerTriangularColumnMatrix::RowHeightPair;
using RowKey = LowerTriangularColumnMatrix::RowKey;

const size_t LowerTriangularColumnMatrix::width() const {
    return width_; 
}

const size_t LowerTriangularColumnMatrix::height() const {
    assert(newMaxHeight_ == maxHeight_);
    return newMaxHeight_ + 1;   // This height should include last row
}

const Key LowerTriangularColumnMatrix::key() const {
   return key_;
}

size_t LowerTriangularColumnMatrix::matrixHeight() const {
    return maxHeight_ + 1;      // This height should include last row
}

// if underlying matrix height is the same as expected matrix height
bool LowerTriangularColumnMatrix::allocated() const {
    return maxHeight_ == newMaxHeight_;
}

// Diagonal flag to indicate if diagonal block should be initialized
// All column matrices should have a diagonal block corresponding to the column index
// In the Cholesky matrix, we want the diagonal block to be on top
// In the Hessian matrix, we don't care about the ordering of blocks
// So we should always keep diagonal block on top
LowerTriangularColumnMatrix::LowerTriangularColumnMatrix(
        const Key key_in, const size_t width_in) 
        : key_(key_in), width_(width_in) {
    // Set up final row first
    preallocateLastRow();

    // Allocate diagonal block
    preallocateBlock(key_, width_, true);
}

bool LowerTriangularColumnMatrix::preallocateBlock(const Key otherKey,
                                               const size_t height,
                                               const bool initialize) {
    auto iterPair = blockStartMap_.insert({otherKey, {newMaxHeight_, height}});
    if(iterPair.second) {
        // If inserted, increase max height and set vector 
        blockStartVec_.back() = {otherKey, {newMaxHeight_, height}};
        newMaxHeight_ += height;

        // Make sure last row is still allocated
        blockStartVec_.push_back({LAST_ROW, {newMaxHeight_, 1}});
        blockStartMap_.at(LAST_ROW) = {newMaxHeight_, 1};

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
void LowerTriangularColumnMatrix::resolveAllocate() {
    assert(newMaxHeight_ >= maxHeight_);    
    if(newMaxHeight_ == maxHeight_) {   // This happens if we don't have anything to allocate
        return;
    }
    if(maxHeight_ == 0) {   
        matrix_ = ColMajorMatrix(newMaxHeight_ + 1, width_);
        matrix_.setZero();
    }
    else {
        // Matrix is not empty, resize without changing values
        matrix_.conservativeResize(newMaxHeight_ + 1, Eigen::NoChange_t());

        // Copy over the last row
        matrix_.block(newMaxHeight_, 0, 1, width_) = matrix_.block(maxHeight_, 0, 1, width_);

        // Set new blocks to 0
        matrix_.block(maxHeight_, 0, newMaxHeight_ - maxHeight_, width_).setZero();
    }
    maxHeight_ = newMaxHeight_;
}

void LowerTriangularColumnMatrix::setZero() {
    matrix_.setZero();
}

void LowerTriangularColumnMatrix::setZero(const Key i) {
    block(i).setZero();
}

// reset blockStart* assignments except for the diagonal block, 
// but don't touch the underlying matrix in case we need that
// memory later. Excess memory will be freed by resize
void LowerTriangularColumnMatrix::resetBlocks() {
    maxHeight_ = 0;
    newMaxHeight_ = 0;
    blockStartVec_.clear();
    blockStartMap_.clear();

    // Set up final row first
    preallocateLastRow();

    // Allocate diagonal block
    preallocateBlock(key_, width_, true);
}

void LowerTriangularColumnMatrix::preallocateLastRow() {
    assert(newMaxHeight_ == 0);
    blockStartMap_.insert({LAST_ROW, {0, 1}});    
    blockStartVec_.push_back({LAST_ROW, {0, 1}});    
}

bool LowerTriangularColumnMatrix::blockExists(const RowKey i) const {
    return blockStartMap_.find(i) != blockStartMap_.end();
}

// Access functions
Block LowerTriangularColumnMatrix::block(const RowKey i) {
    auto pair = blockStartMap_.at(i); 
    return matrix_.block(pair.first, 0, pair.second, width_);
}

const constBlock LowerTriangularColumnMatrix::block(const RowKey i) const {
    auto pair = blockStartMap_.at(i); 
    return matrix_.block(pair.first, 0, pair.second, width_);
}

// Access the underlying matrix with variable height but fixed width
// This is used for when we want to compute the contribution blocks
Block LowerTriangularColumnMatrix::blockRange(
        const size_t startRow, const size_t height) {
    return matrix_.block(startRow, 0, height, width_);
}

const constBlock LowerTriangularColumnMatrix::blockRange(
        const size_t startRow, const size_t height) const {
    return matrix_.block(startRow, 0, height, width_);
}

// Block LowerTriangularColumnMatrix::submatrix(const size_t startRow, const size_t startCol,
//                                          const size_t height, const size_t width) {
//     return matrix_.block(startRow, startCol, height, width);
// }

Block LowerTriangularColumnMatrix::diagonalBlock() {
    assert(maxHeight_ >= width_);
    return matrix_.block(0, 0, width_, width_);
}
const constBlock LowerTriangularColumnMatrix::diagonalBlock() const {
    assert(maxHeight_ >= width_);
    return matrix_.block(0, 0, width_, width_);
}

bool LowerTriangularColumnMatrix::hasBelowDiagonalBlocks() const {
    return newMaxHeight_ > width_; // if matrix doesn't just have the diagonal blocks
}

Block LowerTriangularColumnMatrix::belowDiagonalBlocks() {
    return matrix_.block(width_, 0, maxHeight_ - width_, width_);
}
const constBlock LowerTriangularColumnMatrix::belowDiagonalBlocks() const {
    return matrix_.block(width_, 0, maxHeight_ - width_, width_);
}

const vector<pair<RowKey, RowHeightPair>>& LowerTriangularColumnMatrix::blockStartVec() const {
    return blockStartVec_;
}
const unordered_map<RowKey, RowHeightPair>& LowerTriangularColumnMatrix::blockStartMap() const {
    return blockStartMap_;
}

void LowerTriangularColumnMatrix::reorderBlocks(const vector<Key>& reorderedKeys,
                                                size_t oldLowestKeyIndex) {
    // We can guarantee that the last row will not be changed
    // also guarantee that the order of diagonal block cannot be changed
    assert(reorderedKeys.front() != key_);
    const Key oldLowestKey = blockStartVec_[oldLowestKeyIndex].first;
    size_t firstRow = blockStartMap_[oldLowestKey].first;

    if(2 * firstRow <= newMaxHeight_) {
        // Just replace underlying matrix
        ColMajorMatrix newMatrix(newMaxHeight_ + 1, width_);
        newMatrix.block(0, 0, firstRow, width_) = matrix_.block(0, 0, firstRow, width_);

        size_t curRow = firstRow;
        for(const Key key : reorderedKeys) {
            auto& p = blockStartMap_[key];
            const size_t oldRow = p.first;
            const size_t height = p.second;

            newMatrix.block(curRow, 0, height, width_) = matrix_.block(oldRow, 0, height, width_); 

            p.first = curRow;
            blockStartVec_[oldLowestKeyIndex] = {key, p};
            oldLowestKeyIndex++;

            curRow += height;
        }

        // Assign last row
        newMatrix.block(newMaxHeight_, 0, 1, width_) 
            = matrix_.block(newMaxHeight_, 0, 1, width_); 

        matrix_ = std::move(newMatrix);
    }
    else {
        // Directly change underlying matrix
        ColMajorMatrix newMatrix(newMaxHeight_ - firstRow, width_);

        size_t curRow = 0;
        for(const Key key : reorderedKeys) {
            auto& p = blockStartMap_[key];
            const size_t oldRow = p.first;
            const size_t height = p.second;

            newMatrix.block(curRow, 0, height, width_) = matrix_.block(oldRow, 0, height, width_); 
            p.first = firstRow + curRow;
            blockStartVec_[oldLowestKeyIndex] = {key, p};
            oldLowestKeyIndex++;

            curRow += height;
        }

        matrix_.block(firstRow, 0, newMaxHeight_ - firstRow, width_) = newMatrix;
    }

}


void LowerTriangularColumnMatrix::print(ostream& os) const {
    os << "Column " << key_ << endl << blockStartVec_.size() << endl;
    for(const auto p : blockStartVec_) {
        int rowkey = p.first == LAST_ROW? -1 : p.first;
        size_t row = p.second.first;
        size_t height = p.second.second;
        os << rowkey << " " << row << " " << height << endl;
        os << blockRange(row, height) << "\n\n";
    }
}

void LowerTriangularColumnMatrix::checkInvariant() const {
    assert(blockStartVec_.size() == blockStartMap_.size());
    for(auto p : blockStartVec_) {
        RowKey k = p.first;
        assert(blockStartMap_.find(k) != blockStartMap_.end());
        assert(blockStartMap_.at(k) == p.second);
    }
}

} // namespace gtsam
