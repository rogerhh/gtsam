/**
* @file    JacobianMatrix.cpp
* @brief   Jacobian Matrix
* @author  Roger Hsiao (rogerhh)
* @date    Feb. 8, 2023
*/


#include <gtsam/linear/JacobianMatrix.h>
#include <algorithm>

using namespace std;

namespace gtsam {

using Block = JacobianMatrix::Block;
using constBlock = JacobianMatrix::constBlock;
    
JacobianMatrix::JacobianMatrix() {
    // Allocate b column
    // We're going to use index 0 for the b column
    addColumn(-1, 1);
}

void JacobianMatrix::addColumn(const Key key, const size_t width) {
    // Column indexing starts from 1, because we reserve 0 for column -1
    assert(key + 1 = columns_.size());
    columns_.emplace_back(SparseColumnBlockMatrix(key, width));
}

// FactorIndex
void JacobianMatrix::preallocateBlock(const FactorIndex i, 
                                      const Key key, 
                                      const size_t height, 
                                      bool initialize) {
    assert(key + 1 < columns_.size());
    columns_[key + 1].preallocateBlock(i, height, true);
}

void JacobianMatrix::resolveAllocate(const Key key) {
    assert(key + 1 < columns_.size());
    columns_[key + 1].resolveAllocate();
}

void JacobianMatrix::setZero(const Key key) {
    assert(key + 1 < columns_.size());
    columns_[key + 1].setZero();
}

void JacobianMatrix::setZero(const FactorIndex i, const Key key) {
    assert(key + 1 < columns_.size());
    columns_[key + 1].setZero(i);
}

void JacobianMatrix::resetColumn(const Key key) {
    assert(key + 1 < columns_.size());
    columns_[key + 1].resetBlocks();
}

bool JacobianMatrix::blockExists(const FactorIndex i, const Key key) const {
    return columns_[key + 1].blockExists(i);
}

// Access functions
Block JacobianMatrix::block(const FactorIndex i, const Key key) {
    return columns_[key + 1].block(i);
}
// const JacobianMatrix::constBlock block(const Key i) const {}

// Access the underlying matrix with variable height
// This is used for when we want to compute the contribution blocks
// Block JacobianMatrix::blockRange(const size_t startRow, const size_t height) {}
// const JacobianMatrix::constBlock blockRange(const size_t startRow, const size_t height) const {}

void JacobianMatrix::print(std::ostream& os) {
    os << "Matrix" << endl << columns_.size() << endl;
    for(const auto& col : columns_) {
        col.print(os);
    }
}


}
