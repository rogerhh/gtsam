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

void JacobianMatrix::addColumn(const Key key, const size_t width) {
    columns_.emplace_back(SparseColumnBlockMatrix(key, width, false));
}

// FactorIndex
void JacobianMatrix::preallocateBlock(const FactorIndex i, 
                                      const Key j, 
                                      const size_t height, 
                                      bool initialize) {
    assert(columns_.size() > j);
    columns_[j].preallocateBlock(i, height, true);
}

void JacobianMatrix::resolveAllocate(const Key i) {
    columns_[i].resolveAllocate();
}

void JacobianMatrix::setZero(const Key i) {
    columns_[i].setZero();
}

void JacobianMatrix::setZero(const Key i, const Key j) {
    columns_[j].setZero(i);
}

void JacobianMatrix::resetColumn(const Key i) {
    columns_[i].resetBlocks();
}

bool JacobianMatrix::blockExists(const Key i, const Key j) const {
    return columns_[j].blockExists(i);
}

// Access functions
Block JacobianMatrix::block(const Key i, const Key j) {
    return columns_[j].block(i);
}
// const JacobianMatrix::constBlock block(const Key i) const {}

// Access the underlying matrix with variable height
// This is used for when we want to compute the contribution blocks
// Block JacobianMatrix::blockRange(const size_t startRow, const size_t height) {}
// const JacobianMatrix::constBlock blockRange(const size_t startRow, const size_t height) const {}

void JacobianMatrix::print(std::ostream& os) {
    for(const auto& col : columns_) {
        os << "Column: " << col.key() << endl;
        col.print(os);
    }
}


}
