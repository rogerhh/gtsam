/**
* @file    SparseLowerTriangularBlockMatrix.cpp
* @brief   Access to matrices via blocks of pre-defined sized. Used as a unified Cholesky factor matrix 
* @author  Roger Hsiao (rogerhh)
* @date    Feb. 8, 2023
*/

#include <cassert>
#include <stdexcept>
#include <gtsam/base/SparseLowerTriangularBlockMatrix.h>

using namespace std;

namespace gtsam {

using Block = SparseLowerTriangularBlockMatrix::Block;
using constBlock = SparseLowerTriangularBlockMatrix::constBlock;
using RowHeightPair = SparseLowerTriangularBlockMatrix::RowHeightPair;

void SparseLowerTriangularBlockMatrix::addColumn(const Key key, const size_t width) {
    // For now, you're only allowed to add keys in order
    assert(key == columns_.size());
    if(key != columns_.size()) {
        throw runtime_error("Keys are not added in order!");
    }
    columns_.push_back(LowerTriangularColumnMatrix(key, width));
    needAlloc_.insert(key);
}

size_t SparseLowerTriangularBlockMatrix::colWidth(const Key key) const {
    return columns_[key].width();
}

size_t SparseLowerTriangularBlockMatrix::colHeight(const Key key) const {
    return columns_[key].height();
}

bool SparseLowerTriangularBlockMatrix::preallocateBlock(const Key i, const Key j, bool initialize) {
    if(columns_[j].preallocateBlock(i, colWidth(i), initialize)) {
        needAlloc_.insert(j);
        return true;
    }
    return false;
}

void SparseLowerTriangularBlockMatrix::resolveAllocate() {
    for(const size_t i : needAlloc_) {
        columns_[i].resolveAllocate();
    }
}

void SparseLowerTriangularBlockMatrix::resolveAllocate(const Key key) {
    columns_[key].resolveAllocate();
}

void SparseLowerTriangularBlockMatrix::setZero(const Key i, const Key j) {
    columns_[j].setZero(i);
}

void SparseLowerTriangularBlockMatrix::setZero(const Key i) {
    columns_[i].setZero();
}

void SparseLowerTriangularBlockMatrix::resetColumn(const Key i) {
    columns_[i].resetBlocks();
}

bool SparseLowerTriangularBlockMatrix::blockExists(const Key i, const Key j) const {
    return columns_[j].blockExists(i);
}

Block SparseLowerTriangularBlockMatrix::block(const Key i, const Key j) {
    return columns_[j].block(i);
}

const constBlock SparseLowerTriangularBlockMatrix::block(const Key i, const Key j) const {
    return columns_[j].block(i);
}

Block SparseLowerTriangularBlockMatrix::colBlockRange(
        const Key key, const size_t row, const size_t height) {
    return columns_[key].blockRange(row, height);
}

const constBlock SparseLowerTriangularBlockMatrix::colBlockRange(
        const Key key, const size_t row, const size_t height) const {
    return columns_[key].blockRange(row, height);
}

Block SparseLowerTriangularBlockMatrix::colDiagonalBlock(const Key key) {
    return columns_[key].diagonalBlock();
}

const constBlock SparseLowerTriangularBlockMatrix::colDiagonalBlock(const Key key) const {
    return columns_[key].diagonalBlock();
}

Block SparseLowerTriangularBlockMatrix::colBelowDiagonalBlocks(const Key key) {
    const size_t width = columns_[key].width();
    const size_t height = columns_[key].height();
    return columns_[key].blockRange(width, height - width);
}

const constBlock SparseLowerTriangularBlockMatrix::colBelowDiagonalBlocks(const Key key) const {
    const size_t width = columns_[key].width();
    const size_t height = columns_[key].height();
    return columns_[key].blockRange(width, height - width);
}

LowerTriangularColumnMatrix& SparseLowerTriangularBlockMatrix::column(const Key key) {
    return columns_[key];
}

void SparseLowerTriangularBlockMatrix::checkInvariant() const {

    for(const auto& col : columns_) {
        col.checkInvariant();
        for(const auto p : col.blockStartVec()) {
            const Key k = p.first;
            assert(blockExists(k, col.key()));
            if(k != col.key() && k != LAST_ROW) {
                assert(!blockExists(col.key(), k));
            }
        }
    }
}

void SparseLowerTriangularBlockMatrix::print(std::ostream& os) const {
    os << "Matrix" << endl << columns_.size() << endl;
    for(const auto& col : columns_) {
        col.print(os);
    }
}

void SparseLowerTriangularBlockMatrix::printColumn(std::ostream& os, const Key i) const {
    columns_[i].print(os);
}

} // namespace gtsam

