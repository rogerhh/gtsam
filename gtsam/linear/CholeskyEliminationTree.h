/**
* @file    CholeskyEliminationTree.h
* @brief   Elimination tree structure to perform symbolic factorization and Cholesky factorization
* @author  Roger Hsiao (rogerhh)
* @date    Feb. 8, 2023
*/

#pragma once

#include <gtsam/base/Matrix.h>
#include <gtsam/base/Workspace.h>
#include <gtsam/linear/CholeskyEliminationTreeTypes.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/ISAM2UpdateParams.h>
#include <vector>
#include <unordered_map>
#include <utility>
#include <iostream>
#include <boost/optional.hpp>

namespace gtsam {

class GTSAM_EXPORT CholeskyEliminationTree {
private:
  class Node;
  class Clique;
  class FactorWrapper;   // A wrapper around Factor, can contain nonlinear or linear factors

public:
  typedef std::shared_ptr<Node> sharedNode;
  typedef std::shared_ptr<Clique> sharedClique;
  typedef std::weak_ptr<Clique> weakClique;
  typedef std::shared_ptr<FactorWrapper> sharedFactorWrapper;
  typedef NonlinearFactor::shared_ptr sharedFactor;
  typedef GaussianFactor::shared_ptr sharedLinearFactor;

  static constexpr size_t BLOCK_INDEX_KEY = 0;
  static constexpr size_t BLOCK_INDEX_ROW = 1;
  static constexpr size_t BLOCK_INDEX_HEIGHT = 2;
  
  friend std::ostream& operator<<(std::ostream& os, const Clique& clique);
  friend std::ostream& operator<<(std::ostream& os, const FactorWrapper& factorWrapper);

private:
    std::vector<sharedNode> nodes_;
    // The clique that owns each node. May be shared between nodes. We need this to 
    // avoid circular references and also make sure something owns the sharedClique's
    std::vector<sharedClique> cliques_;   
    sharedClique root_;

    // Transform random key indices into ordered key indices starting at 1
    // 0 is used for the Atb row
    std::unordered_map<Key, RemappedKey> keyTransformMap_;
    // unmappedKeys transforms the sequential remappedKeys back to unmapped keys
    std::vector<Key> unmappedKeys_;

    std::vector<RemappedKey> keyToOrdering_; // Key to index
    std::vector<size_t> orderingToKey_;         // Index to key

    size_t orderingVersion_ = 0;

    struct OrderingLess {
        CholeskyEliminationTree* etree = nullptr;
        OrderingLess() : etree(nullptr) {}
        OrderingLess(CholeskyEliminationTree* etree_in) : etree(etree_in) {}
        bool operator()(const RemappedKey lhs, const RemappedKey rhs) const {
          assert(etree);
          return etree->keyToOrdering_.at(lhs) < etree->keyToOrdering_.at(rhs);
        }
    };

    OrderingLess orderingLess_;

    std::vector<sharedFactorWrapper> factors_;
    std::vector<size_t> factorIndexTransformMap_;

    Workspace workspace_;

    bool postOrder_ = false;
    size_t totalDeltaDim_ = 0;
    
    // TODO: Remove DEBUG
    std::vector<std::vector<RemappedKey>> descendants_;
    std::vector<std::vector<RemappedKey>> changedDescendants_;

public:
  CholeskyEliminationTree();

  void addVariables(const Values& newTheta);
  
  // Mark directly changed keys and keys that we explicitly want to update (extraKeys)
  // observedKeys are the keys associated with the new factors
  void markAffectedKeys(const NonlinearFactorGraph& nonlinearFactors,
                        const FactorIndices& newFactorIndices,
                        const KeySet& relinKeys, 
                        const ISAM2UpdateParams& updateParams,
                        // const boost::optional<FastList<Key>>& extraKeys,
                        RemappedKeySet* affectedKeys,
                        RemappedKeySet* observedKeys);

  // Mark all ancestors of directly changed keys and disconnect child from parent 
  // As parent might change. But do not reset lambdaStructure as that does not change
  void markAncestors(const RemappedKeySet& affectedKeys, RemappedKeySet* markedKeys);

  void updateOrdering(const RemappedKeySet& markedKeys, 
                      const RemappedKeySet& observedKeys,
                      const ISAM2UpdateParams& params);

  void symbolicElimination(const RemappedKeySet& markedKeys);

  void choleskyElimination(const Values& theta);

  void backsolve(VectorValues* delta_ptr, double tol);

  // Reset after each iteration
  void reset();

  // Add all ancestors of unmappedKey to additionalKeys. This is used in marginalization
  void getAffectedKeys(Key unmappedKey, std::set<Key>& additionalKeys) const;

  // Add all descendants of unmappedKey to additionalKeys. This is used in marginalization
  void getAffectedDescendantKeys(Key unmappedKey, std::set<Key>& additionalKeys) const;

  void marginalizeLeaves(
      const FastList<Key>& leafKeys,
      boost::optional<FactorIndices&> marginalFactorsIndices = boost::none,
      boost::optional<FactorIndices&> deletedFactorsIndices = boost::none);
  

  void printOrderingUnmapped(std::ostream& os) const;
  void printOrderingRemapped(std::ostream& os) const;

  // Number of all factors, including deleted factors, marginal factors, and linear factors
  size_t numFactors() const { return factors_.size(); }

  sharedFactor nonlinearFactorAt(FactorIndex factorIndex);

private:
  // Add new unmapped Key to transform map and return the mapped key, 
  // If key already exists, just return the mapped key
  // The new key is added to the end of the variable ordering (but before the last row)
  RemappedKey addRemapKey(const Key unmappedKey);
  RemappedKey getRemapKey(const Key unmappedKey) const;
  Key unmapKey(const RemappedKey remappedKey) const;

  void addNewNode(const Key unmappedKey, const size_t width);
  
  void markKey(const RemappedKey key, RemappedKeySet* markedKeys);
  
  // Perform symbolic elimination on the columns of variable key
  void symbolicEliminateKey(const RemappedKey key);

  void remapConstrainedKeys(const boost::optional<FastMap<Key, int>>& unmappedConstrainedKeys,
                            std::map<RemappedKey, int>* constrainedKeys);

  // Construct CSC matrix of all the nodes denoted by markedKeys
  // markedKeysIndex[key] gives the index of key within all the marked keys
  void constructCSCMatrix(const std::unordered_map<RemappedKey, int>& markedKeysIndex,
                          const std::vector<RemappedKey>& markedKeysVector,
                          int* nEntries,
                          int* nVars,
                          int* nFactors,
                          std::vector<int>* A,
                          std::vector<int>* p);

  void getPartialReordering(const RemappedKeySet& markedKeys, 
                            const RemappedKeySet& observedKeys,
                            const ISAM2UpdateParams& params,
                            std::vector<RemappedKey>* partialOrdering);

  // Post reordering is used to make sure nodes in the same clique is 
  // adjacent in ordering
  // We need this so that when performing a structure altering operation (removing factors
  // or adding factors in the middle, we can retain some structure. Need to prove this)
  void postReordering(const std::vector<RemappedKey>& partialOrdering);
  
  // Allocate required memory for outer products
  // Required memory is the max sum of matrix sizes for any path from root to leaf
  // This function also sets the update mode of the clique to EDIT or RECONSTRUCT
  // We want to minimize clique traversal, so we will set everything whose changed 
  // descendants have a fixed clique to be EDIT, otherwise RECONSTRUCT
  void allocateStack();
  void allocateStackRegular();
  void allocateStackPostOrdering();
  void deallocateStack();

  // Allocate and gather columns of the clique to a new workspace
  // If reconstruct is true, will not gather but simply reset the columns
  void allocateAndGatherClique(sharedClique clique, bool allocate, bool reconstruct);

  void editOrReconstructFromClique(
         sharedClique clique,
         const std::vector<Key>& destCols,
         double sign);
  
  void restoreClique(sharedClique clique);
  
  // Reset the first diagWidth columns in the clique's workspace
  void resetCliqueColumns(sharedClique clique);
  
  // ADD AtA of factors. Handle edits of AtA here as well
  void constructLambdaClique(sharedClique clique, const Values& theta);

  void eliminateClique(sharedClique clique);

  void mergeWorkspaceClique(sharedClique clique);

  // Copy the data in the workspace to a local memory owned by the clique
  void scatterClique(sharedClique clique);

  void backsolveClique(sharedClique clique, 
                       VectorValues* delta_ptr, 
                       double tol);

  bool valuesChanged(const Vector& diff, double tol);

  size_t colWidth(const RemappedKey key) const;

  void marginalizeClique(
      sharedClique clique,
      FactorIndices* deletedFactors,
      FactorIndices* marginalFactors);

  bool checkSortedNoDuplicates(const std::vector<Key>& v) const;

  void checkInvariant_afterSymbolic() const;
  void checkInvariant_afterCholesky() const;
  void checkInvariant_afterBackSolve() const;
  void checkInvariant_afterMarginalize() const;

  // Public testing functions
public:
  void testFactorWrapper(const Values& theta);
};

} // namespace gtsam
