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

  sharedClique parent = nullptr;
  std::set<sharedClique> children;

  // bool hasReconstructAncestor = false;
  bool marked_ = false; 
  // bool backsolve = false;

  MarkedStatus status = UNMARKED;

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

  // Find new parent clique as the lowest nonzero index in any column of the clique
  // not including the diagonal
  void findParent();

  // Set parent to nullptr
  void detachParent();

  void setParent(sharedClique newParent);

  // Merge otherClique into this clique
  void mergeClique(sharedClique otherClique);

  void mergeGatherSources(const std::vector<LocalCliqueColumns>& otherGatherSources);

  // Merge our column structure based on blockIndices into parentColStructure
  void mergeColStructure(std::vector<RemappedKey>* parentColStructure);

  // Populate blockIndices given a colStructure
  void populateBlockIndices(const std::vector<RemappedKey>& colStructure);

  // Reorder underlying matrix. Return false if nothing needs to be done
  bool reorderColumn(BlockIndexVector& newBlockIndices);

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

  void setBacksolve(bool backsolve);

  bool needsBacksolve() const;

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

  void printClique(std::ostream& os);

  // Check if we own our columns
  bool ownsColumns() const;
};

} // namespace gtsam
