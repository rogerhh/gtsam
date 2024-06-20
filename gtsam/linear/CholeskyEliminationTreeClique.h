/**
* @file    CholeskyEliminationTreeClique.h
* @brief   A group of fully connected node in the CholeskyEliminationTree
* @author  Roger Hsiao (rogerhh)
* @date    Feb. 8, 2023
*/

#pragma once

#include <cstddef>
#include <gtsam/linear/CliqueColumns.h>
#include <gtsam/linear/CholeskyEliminationTree.h>
#include <gtsam/linear/CholeskyEliminationTreeTypes.h>
#include <utility>
#include <vector>
#include <unordered_set>
#include <iostream>

namespace gtsam {

class CholeskyEliminationTree::Clique : 
    public std::enable_shared_from_this<CholeskyEliminationTree::Clique> {

private:
  CholeskyEliminationTree* etree = nullptr;

public:

  sharedClique get_ptr() {
    return shared_from_this();
  }

  // Ordered list of nodes belonging to this clique
  std::vector<sharedNode> nodes;

  weakClique parent_;
  std::set<sharedClique> children;

  bool marked_ = false; 
  // bool backsolve = false;

  MarkedStatus status_ = UNMARKED;

  // The memory allocated for this clique
  size_t workspaceIndex = -1;

  // The column structure of this clique, starting from the first column
  BlockIndexVector blockIndices;

  // The ordering version of this clique, if not the same as the etree's ordering version
  // then we need to reorder this clique
  size_t orderingVersion = 0;

  // The clique splitting issue. Currently, if a node in the middle of a clique is marked
  // We would split it into an unmarked clique and marked cliques
  // However, if the unmarked clique is later affected, and it can be merged back into 
  // the marked cliques, we would need to gather column from different locations
  // 06/20/2023: This is needed for deleting factors 
  // or adding factors in the middle of the clique anyways

  // col_data is where the matrix lives during and after cholesky
  // But after markAncestors, cliques might be broken up, col_data_ptr points to the 
  // underlying vector the matrix currently lives in. If col_data_start == 0,
  // then this clique owns the matrix
  // gatherSources is a vector of cliqueColumns, which is where the underlying matrix actually
  // lives. If the size of gatherSources > 1, then the clique's columns can come from
  // multiple previously incontiuous columns, which can be scattered into
  // a contiguous memory though addCliqueColumns()
  std::vector<LocalCliqueColumns> gatherSources;

  // FIXME: make this into data structure in tree traversal
  // How much memory is needed at this clique and all its children
  size_t accumSize = 0;
  // How much memory is needed just for this node
  size_t selfSize = 0;

  Clique(CholeskyEliminationTree* etree);

  // add node to clique
  void addNode(sharedNode node);

  // Mark clique starting from lowest key. Detach all nodes from this clique
  // and add into their own cliques. Children cliques of this clique
  // are kept in this clique. Detach parent and return it
  sharedClique markClique(const RemappedKey lowestKey, RemappedKeySet* markedKeys);

  // Reorder the clique's blocks. Assume the clique's columns are in the same source
  // and we own that source
  void reorderClique();

  sharedClique parent();

  // Find new parent clique as the lowest nonzero index in any column of the clique
  // not including the diagonal
  void findParent();

  // Find new parent after variable reordering
  void reorderAndFindParent();

  // Set parent to nullptr
  void detachParent();

  void setParent(sharedClique newParent);

  // Set parent without detaching the original parent
  // Assuming that the original parent will be destroyed
  void setParent2(sharedClique newParent);

  // Merge parent otherClique into this clique
  void mergeClique(sharedClique otherClique);

  void mergeGatherSources(const std::vector<LocalCliqueColumns>& otherGatherSources);
  
  // Merge childClique into this clique
  // This is to support merging multiple branches
  void mergeClique2(sharedClique childClique);

  void mergeGatherSources2(std::vector<LocalCliqueColumns>& childGatherSources);

  // Merge our column structure based on blockIndices into parentColStructure
  void mergeColStructure(std::vector<RemappedKey>* parentColStructure);

  // Populate blockIndices given a colStructure
  void populateBlockIndices(const std::vector<RemappedKey>& colStructure);

  // Extend a blockIndices with a colStructure
  void extendBlockIndices(const std::vector<RemappedKey>& colStructure);

  // Returns true if clique contributes to an ancestor that is marked
  // Only need to check highest key that is not 0
  bool hasMarkedAncestor();

  // Determine if a clique needs to be edited or reconstructed
  // Right now, the algorithm should be that, if any changed descendant is unmarked
  // we set it to EDIT
  void setEditOrReconstruct();

  void checkEditOrReconstruct(
      MarkedStatus mode, std::vector<Key>* destCols);

  // After splitting a clique, if there is a reordering, we need to deallocate 
  // all old columns as all the marked variables will be reconstruct
  void deallocateExcessCols();

  // Iterate over all nodes and set their MarkedStatus
  // If a node is NEW and the clique is EDIT, the node remains NEW
  void setNodeStatus();
  void setNodeStatusMarginalize();

  const MarkedStatus& status() const { return status_; }

  void setStatusUnmarked() {
    assert(status_ == UNMARKED || status_ == EDIT || status_ == RECONSTRUCT);
    status_ = UNMARKED;
  }

  void setStatusEdit() {
    assert(status_ == UNMARKED || status_ == EDIT || status_ == RECONSTRUCT);
    status_ = EDIT;
  }

  void setStatusReconstruct() {
    assert(status_ == UNMARKED || status_ == RECONSTRUCT);
    status_ = RECONSTRUCT;
  }
  void setStatusMarginalized() {
    assert(status_ == UNMARKED);
    status_ = MARGINALIZED;
  }

  void setBacksolve(bool backsolve);

  bool needsBacksolve() const;

  // Split a clique into 2, where the second clique owns its own data. This is used 
  // in marginalization
  sharedClique splitClique(size_t splitIndex);

  // Reset all member variables after elimination and all nodes
  void resetAfterCholesky();

  void resetAfterBacksolve();

  sharedNode front();
  sharedNode back();

  RemappedKey frontKey() const;
  RemappedKey backKey() const;

  // Return true if this clique is associated with the last row (key 0)
  bool isLastRow() const;

  bool marked() const;

  // Number of nodes in this clique
  size_t cliqueSize() const;

  size_t diagonalHeight() const;

  size_t subdiagonalHeight() const;
  //
  // Height of the matrix 
  size_t height() const;

  // Width of the matrix
  size_t width() const;

  void printBlockIndices(std::ostream& os);

  void printClique(std::ostream& os);

  // Check if we own our columns
  bool ownsColumns() const;

  void checkInvariant() const;

  bool match(const std::vector<RemappedKey>& v) const {
    if(v.size() != blockIndices.size()) {
      return false;
    }
    for(int i = 0; i < v.size(); i++) {
        auto& [key, row, height] = blockIndices[i];
        std::cout << key << std::endl;
        if(v[i] != key) {
            return false;
        }
    }
    return true;
  }

  void deleteClique();

  /**************
   * Cost computation
   *************/
  int64_t SYM_COST = 2000;
  int64_t SYM_COST_REORDER = 4000;

  bool parallelizable = true;
  CostStatus costStatus = COST_UNMARKED;
  CostStatus nextCostStatus = COST_UNMARKED;
  int64_t markedCost = -1;
  int64_t fixedCost = -1;
  int64_t restoreCost = -1;
  int64_t symCost = -1;
  int64_t AtACost = -1;
  int64_t syrkCost = -1;
  int64_t addCost = -1;
  int64_t cholCost = -1;
  int64_t backsolveCost = -1;

  void resetCost();
  // Even though each clique is assumed to run under 1 thread
  // We still split the cost of the clique by num_threads except for the root clique
  int64_t computeCostMarked(int num_threads);
  int64_t computeCostFixed(int num_threads);
  int64_t computeCostBacksolve(int num_threads);

  int memSize = 0;
  int subdiagSize = 0;
  int factorSize = 0;
  int lastRelinStep = 0;

public:
  static uint64_t merge1_cycles;
  static uint64_t merge2_cycles;
  static uint64_t merge3_cycles;
  static uint64_t merge4_cycles;
  static uint64_t merge5_cycles;
  static uint64_t merge6_cycles;

};

} // namespace gtsam
