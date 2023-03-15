/**
* @file    SparseSymmetricBlockMatrix.cpp
* @brief   Represents a column matrix of fixed width and predefined block heights
* @author  Roger Hsiao (rogerhh)
* @date    Feb. 8, 2023
*/

#include <cassert>
#include <stdexcept>
#include <gtsam/base/SparseSymmetricBlockMatrix.h>

using namespace std;

namespace gtsam {

using Block = SparseSymmetricBlockMatrix::Block;
using constBlock = SparseSymmetricBlockMatrix::constBlock;
using RowHeightPair = SparseSymmetricBlockMatrix::RowHeightPair;

void SparseSymmetricBlockMatrix::addColumn(const Key key, const size_t width) {
    // For now, you're only allowed to add keys in order
    assert(key == columns_.size());
    if(key != columns_.size()) {
        throw runtime_error("Keys are not added in order!");
    }
    columns_.push_back(SparseColumnBlockMatrix(key, width));
    needAlloc_.insert(key);
}

size_t SparseSymmetricBlockMatrix::colWidth(const Key key) const {
    return columns_[key].width();
}

size_t SparseSymmetricBlockMatrix::colHeight(const Key key) const {
    return columns_[key].height();
}

void SparseSymmetricBlockMatrix::preallocateBlock(const Key i, const Key j, bool initialize) {
    if(i >= j) {
        if(columns_[j].preallocateBlock(i, colWidth(i), initialize)) {
            needAlloc_.insert(j);
        }
    }
    else {
        if(columns_[i].preallocateBlock(j, colWidth(j), initialize)) {
            needAlloc_.insert(i);
        }
    }
}

void SparseSymmetricBlockMatrix::resolveAllocate() {
    for(const size_t i : needAlloc_) {
        columns_[i].resolveAllocate();
    }
}

void SparseSymmetricBlockMatrix::setZero(const Key i, const Key j) {
    if(i >= j) {
        columns_[j].setZero(i);
    }
    else {
        columns_[i].setZero(j);
    }
}

void SparseSymmetricBlockMatrix::setZero(const Key i) {
    columns_[i].setZero();
}

void SparseSymmetricBlockMatrix::resetColumn(const Key i) {
    columns_[i].resetBlocks();
}

bool SparseSymmetricBlockMatrix::blockExists(const Key i, const Key j) const {
    if(i >= j) {
        return columns_[j].blockExists(i);
    }
    else {
        return columns_[i].blockExists(j);
    }
}

Block SparseSymmetricBlockMatrix::block(const Key i, const Key j) {
    if(i >= j) {
        return columns_[j].block(i);
    }
    else {
        return columns_[i].block(j);
    }
}

const constBlock SparseSymmetricBlockMatrix::block(const Key i, const Key j) const {
    if(i >= j) {
        return columns_[j].block(i);
    }
    else {
        return columns_[i].block(j);
    }
}

Block SparseSymmetricBlockMatrix::colBlockRange(
        const Key key, const size_t row, const size_t height) {
    return columns_[key].blockRange(row, height);
}

const constBlock SparseSymmetricBlockMatrix::colBlockRange(
        const Key key, const size_t row, const size_t height) const {
    return columns_[key].blockRange(row, height);
}

void SparseSymmetricBlockMatrix::checkInvariant() const {

    for(const auto& col : columns_) {
        col.checkInvariant();
        for(const auto p : col.blockStartVec()) {
            const Key k = p.first;
            assert(col.key() <= k);
        }
    }
}

void SparseSymmetricBlockMatrix::print(std::ostream& os) const {
    for(const auto& col : columns_) {
        os << "Column: " << col.key() << endl;
        col.print(os);
    }
}

void SparseSymmetricBlockMatrix::printColumn(std::ostream& os, const Key i) const {
    columns_[i].print(os);
}

} // namespace gtsam
