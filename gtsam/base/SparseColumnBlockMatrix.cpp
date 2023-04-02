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

SparseColumnBlockMatrix::SparseColumnBlockMatrix(
        const Key key_in, const size_t width_in) 
        : key_(key_in), width_(width_in) {}

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
        matrix_ = ColMajorMatrix(newMaxHeight_, width_);
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

// reset blockStart* assignments 
// but don't touch the underlying matrix in case we need that
// memory later. Excess memory will be freed by resize
void SparseColumnBlockMatrix::resetBlocks() {
    maxHeight_ = 0;
    newMaxHeight_ = 0;
    blockStartVec_.clear();
    blockStartMap_.clear();
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

const vector<pair<Key, RowHeightPair>>& SparseColumnBlockMatrix::blockStartVec() const {
    return blockStartVec_;
}
const unordered_map<Key, RowHeightPair>& SparseColumnBlockMatrix::blockStartMap() const {
    return blockStartMap_;
}

void SparseColumnBlockMatrix::reorderBlocks(const std::vector<RowKey>& reorderedKeys, 
                                            size_t oldLowestKeyIndex) {
    
    // We can guarantee that the last row will not be changed
    const RowKey oldLowestKey = blockStartVec_[oldLowestKeyIndex].first;
    size_t firstRow = blockStartMap_[oldLowestKey].first;

    if(2 * firstRow <= newMaxHeight_) {
        // Just replace underlying matrix
        ColMajorMatrix newMatrix(newMaxHeight_, width_);
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

void SparseColumnBlockMatrix::print(ostream& os) const {
    int intkey = key_ == -1? -1 : key_;
    os << "Column " << intkey << endl << blockStartVec_.size() << endl;
    for(const auto p : blockStartVec_) {
        RowKey rowkey = p.first;
        size_t row = p.second.first;
        size_t height = p.second.second;
        os << rowkey << ": " << blockRange(row, height).transpose() << endl;
        // os << rowkey << " " << row << " " << height << endl;
        // os << blockRange(row, height) << "\n\n";
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
