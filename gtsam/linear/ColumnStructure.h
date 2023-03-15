/**
* @file    ColumnStructure.h
* @brief   Symbolic column structure
* @author  Roger Hsiao (rogerhh)
* @date    Feb. 8, 2023
*/

#pragma once

#include <gtsam/base/Matrix.h>
#include <gtsam/base/MatrixSerialization.h>
#include <gtsam/base/FastVector.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/linear/CholeskyEliminationTree.h>
#include <map>

namespace gtsam {

class CholeskyEliminationTree::ColumnStructure {
    std::map<Key, size_t> col;
};

} // namespace gtsam
