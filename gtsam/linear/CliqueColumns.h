/**
* @file    CliqueColumns.h
* @brief   Representation matrix data in the columns of a clique. There are 2 different types of clique columns, LocalCliqueColumns and WorkspaceCliqueColumns. A LocalCliqueColumns owns its own matrix. A WorkspaceCliqueColumns only has a view on the matrix and does not own its own data, and is just needed to implement addCliqueColumns.
* @author  Roger Hsiao (rogerhh)
* @date    Jun. 25, 2023
*/

#pragma once

#include <cstddef>
#include <gtsam/linear/CholeskyEliminationTreeTypes.h>
#include <utility>
#include <vector>
#include <unordered_set>
#include <iostream>
#include <cassert>
#include <cstring>

extern "C" {
#include <gtsam/linear/gemmini_functions.h>
}

namespace gtsam {

class CliqueColumns {

/*
 *  A CliqueColumns stores the entries of the Cholesky factor. matrix gives a view on the underlying matrix; however, the CliqueColumns doesn't have to own the data. The blockIndices stores the variable keys, ordering, and widths associated with the clique. As the Cholesky factor is symmetric, we only store the subdiagonal part. However, since we might need to split a CliqueColumns into multiple CliqueColumns, the startIndex and endIndex member variables denote the width of the CliqueColumns. 
 * */

public:

  struct Scatter {
    size_t srcRow = 0;
    size_t destRow = 0;
    size_t height = 0;
  };

protected:
  // Pointer to the underlying data of the matrix. Construct an Eigen::Map whenever a view of the matrix is needed
  GEMMINI_TYPE* origMatrixData_ = nullptr;
  GEMMINI_TYPE* matrixData_ = nullptr;     // Pointer to the start of this column
  size_t rows_ = 0, cols_ = 0;

  BlockIndexVector* blockIndices_;

  size_t startIndex_ = 0, endIndex_ = -1;

  inline void vectorAdd(GEMMINI_TYPE* dest, GEMMINI_TYPE* src, size_t len) {
    for(size_t i = 0; i < len; i++) {
      *(dest + i) += *(src + i);
    }
  }

public:

  CliqueColumns() {}

  CliqueColumns(GEMMINI_TYPE* origMatrixData_in, 
                BlockIndexVector* blockIndices_in,
                size_t startIndex_in,
                size_t endIndex_in) 
  : origMatrixData_(origMatrixData_in), 
    blockIndices_(blockIndices_in), 
    startIndex_(startIndex_in), endIndex_(endIndex_in) {
    assert(std::get<BLOCK_INDEX_KEY>(blockIndices_->back()) == 0);
    rows_ = std::get<BLOCK_INDEX_ROW>(blockIndices_->back()) + 1;
    size_t c1 = std::get<BLOCK_INDEX_ROW>(blockIndices_->at(startIndex_));
    size_t c2 = std::get<BLOCK_INDEX_ROW>(blockIndices_->at(endIndex_ - 1)) 
                + std::get<BLOCK_INDEX_HEIGHT>(blockIndices_->at(endIndex_ - 1));
    cols_ = c2 - c1;

    assert(c2 > c1);

    matrixData_ = origMatrixData_ + rows_ * c1;
  }

  CliqueColumns(GEMMINI_TYPE* matrixData_in, 
                BlockIndexVector* blockIndices_in,
                size_t startIndex_in)
  : CliqueColumns(matrixData_in, blockIndices_in, startIndex_in, blockIndices_in->size()) {}


  CliqueColumns(GEMMINI_TYPE* matrixData_in, 
                BlockIndexVector* blockIndices_in) 
  : CliqueColumns(matrixData_in, blockIndices_in, 0, blockIndices_in->size()) {}

  CliqueColumns(const CliqueColumns& other) = default;

  // Getters
  bool allocated() const {
    return origMatrixData_ != nullptr;
  }

  size_t rows() const {
    return rows_;
  }

  size_t cols() const {
    return cols_;
  }

  size_t startIndex() const {
    return startIndex_;
  }

  size_t endIndex() const {
    return endIndex_;
  }

  Eigen::Map<ColMajorMatrix<GEMMINI_TYPE>> matrix() {
    return Eigen::Map<ColMajorMatrix<GEMMINI_TYPE>>(matrixData_, rows_, cols_);
  }

  CliqueColumns split(size_t splitIndex);

  // Merge the other CliqueColumns into this CliqueColumns
  // iff their matrixData is contiguous in memory
  bool merge(const CliqueColumns& other);

  // Get the row mappings from one CliqueColumns to another
  // The blockIndices of the src CliqueColumns must be a subset of this one
  // and must be in the same order
  std::vector<Scatter> getScatterMap(CliqueColumns& src);

  // Get the row mappings from one CliqueColumns to another
  // The blockIndices of the src CliqueColumns must be a subset of this one
  // But they don't have to be in the same order
  std::vector<Scatter> getScatterMapReordered(CliqueColumns& src);

  // Add the entries of the CliqueColumn src into the entries of this. As the src CliqueColumns blockIndices is the subset 
  // of this->blockIndices, a scatter operation is needed.
  void addCliqueColumns(CliqueColumns& src, bool reordered=false);

  // Add the entries of the src CliqueColumns into the entries of this with the help
  // of a row mapping from the src CliqueColumns to this one
  // This is for the possible future extension of spliting the workspace matrix into
  // 2 parts, but only needing 1 scatterMap
  void addCliqueColumns(CliqueColumns& src, const std::vector<Scatter>& scatterMap);

  friend std::ostream& operator<<(std::ostream& os, const CliqueColumns& cliqueColumns);
 
};

class LocalCliqueColumns : public CliqueColumns {

/* Multiple LocalCliqueColumns may shared the same underlying data matrixSource, although the LocalCliqueColumns whose startIndex_=0 is the true owner of the data
 * */

private:
    std::shared_ptr<std::vector<GEMMINI_TYPE>> matrixSource_ = nullptr;
    std::shared_ptr<BlockIndexVector> blockIndicesSource_ = nullptr;

public:
    LocalCliqueColumns(std::shared_ptr<std::vector<GEMMINI_TYPE>> matrixSource_in, 
                       std::shared_ptr<BlockIndexVector> blockIndicesSource_in,
                       size_t startIndex_in, 
                       size_t endIndex_in) 
    : CliqueColumns(matrixSource_in->data(), blockIndicesSource_in.get(), 
                    startIndex_in, endIndex_in),
      matrixSource_(matrixSource_in), blockIndicesSource_(blockIndicesSource_in) {}

    LocalCliqueColumns(std::shared_ptr<std::vector<GEMMINI_TYPE>> matrixSource_in, 
                       std::shared_ptr<BlockIndexVector> blockIndicesSource_in,
                       size_t startIndex_in)
    : LocalCliqueColumns(matrixSource_in, blockIndicesSource_in, 
                         startIndex_in, blockIndicesSource_in->size()) {} 

    LocalCliqueColumns(const LocalCliqueColumns& other) = default;

    LocalCliqueColumns split(size_t splitIndex) {
      LocalCliqueColumns res(matrixSource_, blockIndicesSource_, splitIndex, endIndex_);     
      endIndex_ = splitIndex;
      cols_ = cols_ - res.cols_;
      return res;
    }

    bool ownsData() const {
      return startIndex_ == 0;
    }

    // Dealocate data that we own
    void deallocateSelfData() {
      assert(matrixSource_->size() >= rows_ * cols_);
      matrixSource_->resize(matrixSource_->size() - rows_ * cols_);
    }

    size_t getEndCol() const {
      size_t endCol = std::get<BLOCK_INDEX_ROW>(blockIndicesSource_->at(endIndex_ - 1)) 
                      + std::get<BLOCK_INDEX_HEIGHT>(blockIndicesSource_->at(endIndex_ - 1));
      return endCol;
    }

    bool isEndOfCliqueColumns() const {
      // This clique columns has to be the last clique columns in the matrix source
      return getEndCol() * rows_ == matrixSource_->size();
    }

    void forceOwn() {
      if(ownsData()) { return; }

      // This clique columns has to be the last clique columns in the matrix source
      assert(isEndOfCliqueColumns());

      auto oldMatrixSource = matrixSource_;
      auto oldMatrix = matrix();
      auto oldRows = rows_;
      auto oldCols = cols_;
      auto newMatrixSource = std::make_shared<std::vector<GEMMINI_TYPE>>(rows_ * cols_, 0);
      auto newBlockIndicesSource = std::make_shared<BlockIndexVector>();
      auto itOffset = blockIndicesSource_->begin() + startIndex_;
      const size_t rowOffset = std::get<BLOCK_INDEX_ROW>(*itOffset);
      newBlockIndicesSource->reserve(blockIndicesSource_->size() - startIndex_);
      for(auto it = itOffset; it != blockIndicesSource_->end(); it++) {
        const auto&[key, row, height] = *it;
        newBlockIndicesSource->push_back({key, row - rowOffset, height});
      }
      endIndex_ -= startIndex_;
      startIndex_ = 0;

      *this = LocalCliqueColumns(newMatrixSource, newBlockIndicesSource, startIndex_, endIndex_);
  
      assert(rows_ == oldRows - rowOffset);

      auto newMatrix = matrix();
      Eigen::Block<Eigen::Map<ColMajorMatrix<GEMMINI_TYPE>>> 
        oldBlock(oldMatrix, rowOffset, 0, rows_, cols_);
      Eigen::Block<Eigen::Map<ColMajorMatrix<GEMMINI_TYPE>>> 
        newBlock(newMatrix, 0, 0, rows_, cols_);
      newBlock = oldBlock;

      assert(ownsData());

      // Release old matrix source memtory
      oldMatrixSource->resize(oldMatrixSource->size() - oldCols * oldRows);

    }

    std::shared_ptr<BlockIndexVector> blockIndicesSource() { return blockIndicesSource_; }

    GEMMINI_TYPE* get_ptr() {
      assert(ownsData());

      return matrixSource_->data();
    }
};

}   // namespace gtsam
