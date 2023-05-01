/**
* @file    LowerTriangularColumnMatrix.cpp
* @brief   Represents a column for a lower triangular matrix, contains a diagonal block and a row at the end for Atb. Column blocks need to support insertion in the middle
* @author  Roger Hsiao (rogerhh)
* @date    Mar. 15, 2023
*/

#include <cassert>
#include <unordered_map>
#include <gtsam/base/LowerTriangularColumnMatrix.h>

using namespace std;

namespace gtsam {

using Block = LowerTriangularColumnMatrix::Block;
using constBlock = LowerTriangularColumnMatrix::constBlock;
using RowHeightPair = LowerTriangularColumnMatrix::RowHeightPair;
using RowKey = LowerTriangularColumnMatrix::RowKey;
using BlockIndexVector = LowerTriangularColumnMatrix::BlockIndexVector;

const size_t LowerTriangularColumnMatrix::LAST_ROW = -1;

const size_t LowerTriangularColumnMatrix::width() const {
    return width_; 
}

const size_t LowerTriangularColumnMatrix::height() const {
    return newMaxHeight_ + 1;   // This height should include last row
}

size_t LowerTriangularColumnMatrix::matrixHeight() const {
    // if(maxHeight_ + 1 != matrix_.rows()) {
    //     cout << "Max height wrong! " << maxHeight_ << " " << matrix_.rows() << endl;
    // }
    assert(maxHeight_ + 1 == matrix_.rows());
    return maxHeight_ + 1;      // This height should include last row
}

const Key LowerTriangularColumnMatrix::key() const {
   return key_;
}

// if underlying matrix height is the same as expected matrix height
bool LowerTriangularColumnMatrix::allocated() const {
    return maxHeight_ != 0 && maxHeight_ == newMaxHeight_;
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
}

// void LowerTriangularColumnMatrix::preallocateBlocks(const vector<pair<RowKey, size_t>>& blocks) {
void LowerTriangularColumnMatrix::preallocateBlocks(const BlockIndexVector& blocks,
                                                    size_t startIndex) {
    blockIndices_.resize(blocks.size() - startIndex);

    size_t startRow = blocks[startIndex].second.first;
    for(size_t i = 0; i < blockIndices_.size(); i++) {
        blockIndices_[i] = blocks[i + startIndex];
        blockIndices_[i].second.first -= startRow;
    }

    assert(blockIndices_.back().first == -1);
    newMaxHeight_ = blockIndices_.back().second.first;
}


// void LowerTriangularColumnMatrix::preallocateBlocks(const vector<pair<RowKey, size_t>>& blocks) {
void LowerTriangularColumnMatrix::preallocateBlocks(const vector<pair<RowKey, size_t>>& blocks,
                                                    size_t startIndex) {
    int i = 0;
    if(blockIndices_.size() >= 2) {
        Key lastOldKey = blockIndices_[blockIndices_.size() - 2].first;

        // Find the first old block and add new blocks behind that block
        for(i = blocks.size() - 1; i >= startIndex; i--) {
            // if(blocks[i].first == blockIndices_.back().first) {
            if(blocks[i].first == lastOldKey) {
                break;
            }
        }
        i++;
    }
    size_t index = blockIndices_.size() - 1;    // Add new blocks before the -1 row
    blockIndices_.resize(blockIndices_.size() + blocks.size() - i);
    for(; i < blocks.size(); i++) {
        // cout << index << " " << i << endl;
        blockIndices_[index] = {blocks[i].first, {newMaxHeight_, blocks[i].second}};
        newMaxHeight_ += blocks[i].second;
        index++;
    }
    blockIndices_.back() = {LAST_ROW, {newMaxHeight_, 1}};
    // cout << "newMaxHeight_ = " << newMaxHeight_ << " New block indices: ";
    // for(auto p : blockIndices_) {
    //     cout << "[" << p.first << " " << p.second.first << " " << p.second.second << "] ";
    // }
    // cout << endl;
}

void LowerTriangularColumnMatrix::preallocateBlock(const pair<RowKey, size_t>& block) {
    blockIndices_.push_back({block.first, {newMaxHeight_, block.second}});
    newMaxHeight_ += block.second;
    assert(0);
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
    assert(maxHeight_ + 1 == matrix_.rows());
}

void LowerTriangularColumnMatrix::setZero() {
    matrix_.setZero();
}

// reset blockStart* assignments except for the diagonal block, 
// but don't touch the underlying matrix in case we need that
// memory later. Excess memory will be freed by resize
void LowerTriangularColumnMatrix::resetBlocks() {
    maxHeight_ = 0;
    newMaxHeight_ = 0;
    blockIndices_.clear();

    // Set up final row first
    preallocateLastRow();
}

void LowerTriangularColumnMatrix::preallocateLastRow() {
    assert(newMaxHeight_ == 0);
    blockIndices_.push_back({LAST_ROW, {0, 1}});    
}

Block LowerTriangularColumnMatrix::lastRow() {
    assert(maxHeight_ + 1 == matrix_.rows());
    return matrix_.block(maxHeight_, 0, 1, width_);
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

const BlockIndexVector& LowerTriangularColumnMatrix::blockIndices() const {
    return blockIndices_;
}

void LowerTriangularColumnMatrix::reorderBlocks(const vector<Key>& reorderedKeys,
                                                size_t oldLowestKeyIndex) {
    unordered_map<Key, RowHeightPair> indexMap;
    for(size_t i = oldLowestKeyIndex; i < blockIndices_.size(); i++) {
        indexMap.insert(blockIndices_[i]);
    }

    size_t firstRow = blockIndices_[oldLowestKeyIndex].second.first;
    if(firstRow * 2 >= maxHeight_) {
        // Just replace underlying matrix
        ColMajorMatrix newMatrix(maxHeight_ + 1, width_);
        newMatrix.block(0, 0, firstRow, width_) = matrix_.block(0, 0, firstRow, width_);

        size_t curRow = firstRow;
        for(const Key key : reorderedKeys) {
            auto& p = indexMap.at(key);
            const size_t oldRow = p.first;
            const size_t height = p.second;

            newMatrix.block(curRow, 0, height, width_) = matrix_.block(oldRow, 0, height, width_); 

            p.first = curRow;
            blockIndices_[oldLowestKeyIndex] = {key, p};
            oldLowestKeyIndex++;

            curRow += height;
        }

        // Copy last row
        newMatrix.block(maxHeight_, 0, 1, width_) = matrix_.block(maxHeight_, 0, 1, width_);

        matrix_ = std::move(newMatrix);
    
    }
    else {
        // Directly change underlying matrix
        const size_t remainingHeight = newMaxHeight_ - firstRow;
        ColMajorMatrix newMatrix(remainingHeight, width_);

        size_t curRow = 0;
        for(const Key key : reorderedKeys) {
            auto& p = indexMap.at(key);
            const size_t oldRow = p.first;
            const size_t height = p.second;

            newMatrix.block(curRow, 0, height, width_) = matrix_.block(oldRow, 0, height, width_); 
            p.first = firstRow + curRow;
            blockIndices_[oldLowestKeyIndex] = {key, p};
            oldLowestKeyIndex++;

            curRow += height;

        }

        matrix_.block(firstRow, 0, remainingHeight, width_) = newMatrix;
    }

    assert(maxHeight_ + 1 == matrix_.rows());
}


void LowerTriangularColumnMatrix::print(ostream& os) const {
    os << "Column " << key_ << endl << blockIndices_.size() << endl;
    for(const auto p : blockIndices_) {
        int rowkey = p.first == LAST_ROW? -1 : p.first;
        size_t row = p.second.first;
        size_t height = p.second.second;
        os << rowkey << " " << row << " " << height << endl;
        os << blockRange(row, height) << "\n\n";
    }
}

void LowerTriangularColumnMatrix::printBlockIndices(std::ostream& os) const {
    os << "Column " << key_ << " block indices: ";
    for(const auto p : blockIndices_) {
        os << "[" << p.first << " " << p.second.first << " " << p.second.second << "] ";
    }
    os << endl;
}

void LowerTriangularColumnMatrix::checkInvariant() const {
    unordered_set<RowKey> rows;
    size_t heightSum = 0;
    for(auto p : blockIndices_) {
        auto inserted = rows.insert(p.first);
        assert(inserted.second == true);
        assert(p.second.first == heightSum);
        heightSum += p.second.second;
    }
    assert(blockIndices_.back().first == LAST_ROW);
    assert(heightSum == newMaxHeight_ + 1);

    if(allocated()) {
        assert(newMaxHeight_ == maxHeight_);
    }
}

} // namespace gtsam
