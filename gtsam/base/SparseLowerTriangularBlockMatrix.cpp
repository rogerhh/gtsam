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
    columns_.push_back(SparseColumnBlockMatrix(key, width));
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

SparseColumnBlockMatrix& SparseLowerTriangularBlockMatrix::column(const Key key) {
    return columns_[key];
}

void SparseLowerTriangularBlockMatrix::checkInvariant() const {

    for(const auto& col : columns_) {
        col.checkInvariant();
        for(const auto p : col.blockStartVec()) {
            const Key k = p.first;
            assert(blockExists(k, col.key()));
            if(k != col.key()) {
                assert(!blockExists(col.key(), k));
            }
        }
    }
}

void SparseLowerTriangularBlockMatrix::print(std::ostream& os) const {
    for(const auto& col : columns_) {
        os << "Column: " << col.key() << endl;
        col.print(os);
    }
}

void SparseLowerTriangularBlockMatrix::printColumn(std::ostream& os, const Key i) const {
    columns_[i].print(os);
}

} // namespace gtsam

// #include <gtsam/base/SparseUpperTriangularBlockMatrix.h>
// 
// namespace gtsam {
// 
// using SparseUpperTriangularBlockMatrix::Block;
// using SparseUpperTriangularBlockMatrix::constBlock;
// 
// SparseUpperTriangularBlockMatrix::ColumnMatrix::ColumnMatrix(
//         const Key key_in, const size_t width_in)
//         : key(key_in), width(width_in) {
//     
//     } 
// 
// bool SparseUpperTriangularBlockMatrix::ColumnMatrix::preallocateOrInitialize(
//         const Key otherKey, 
//         const size_t height,
//         const bool initialize) {
//     assert(otherKey <= key);
// 
//     auto iterPair = blockStart.insert({otherKey, {0, 0}});
//     if(iterPair.second) {
//         // if inserted, then requested block did not exist before
//         iterPair.first->second.first = newMaxHeight;
//         iterPair.first->second.second = height;
//         newMaxHeight += height;
//     }
//     else if(initialize) {
//         assert(iterPair.first->second.second == height); 
//         const size_t blockStartRow = iterPair.first->second.first;
//         matrix.block(blockStartRow, 0, height, width).setZero();
//     }
//     return iterPair.second;
// }
// 
// void SparseUpperTriangularBlockMatrix::ColumnMatrix::resolveAllocate() {
//     assert(newMaxHeight > maxHeight);    
//     if(maxHeight == 0) {
//         matrix = RowMajorMatrix(newMaxHeight, width);
//     }
//     else {
//         matrix.conservativeResize(newMaxHeight, Eigen::NoChange_t());
//     }
//     matrix.block(maxHeight, 0, newMaxHeight - maxHeight, width).setZero();
//     maxHeight = newMaxHeight;
// }
// 
// void SparseUpperTriangularBlockMatrix::ColumnMatrix::resetBlocks() {
//     blockStart.clear();
//     maxHeight = 0;
//     newMaxHeight = 0;
// }
// 
// Block SparseUpperTriangularBlockMatrix::ColumnMatrix::block(const Key i) {
//     assert(i >= key);
//     auto pair = blockStart.at(i); 
//     return matrix.block(pair.first, 0, pair.second, width);
// }
// 
// const constBlock SparseUpperTriangularBlockMatrix::ColumnMatrix::block(const Key i) const {
//     assert(i >= key);
//     auto pair = blockStart.at(i); 
//     return matrix.block(pair.first, 0, pair.second, width);
// }
// 
// void SparseUpperTriangularBlockMatrix::addColumn(const Key key, const size_t width) {
//     assert(key == columnMatrices_.size());
//     columnMatrices_.push_back(ColumnMatrix(key, width));
// }
// 
// void SparseUpperTriangularBlockMatrix::preallocateOrInitialize(
//         const Key i, const Key j, const bool initialize) {
//     // assert(i <= j);  // Getting rid of this 
//     if(i <= j) {
//         const size_t height = columnMatrices_[i].width;
//         bool alloc = columnMatrices_[j].preallocateOrInitialize(i, height, initialize);
//         if(alloc) {
//             needAlloc.insert(j);
//         }
//     }
//     else {
//         const size_t height = columnMatrices_[j].width;
//         bool alloc = columnMatrices_[i].preallocateOrInitialize(j, height, initialize);
//         if(alloc) {
//             needAlloc.insert(i);
//         }
//     }
// }
// 
// void SparseUpperTriangularBlockMatrix::resolveAllocate() {
//     for(const size_t i : needAlloc) {
//         columnMatrices_[i].resolveAllocate();
//     }
//     needAlloc.clear();
// }
// 
// void SparseUpperTriangularBlockMatrix::resetColumn(const Key key) {
//     columnMatrices_[key].resetBlocks();
// }
// 
// Block SparseUpperTriangularBlockMatrix::block(const Key i, const Key j) {
//     // assert(i <= j);
//     return (i <= j)? columnMatrices_[j].block(i) : columnMatrices_[i].block(j);
// }
// 
// const constBlock SparseUpperTriangularBlockMatrix::block(const Key i, const Key j) const {
//     // assert(i <= j);
//     return (i <= j)? columnMatrices_[j].block(i) : columnMatrices_[i].block(j);
// }
// 
// } // namespace gtsam
