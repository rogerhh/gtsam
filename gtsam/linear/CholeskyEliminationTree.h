/**
* @file    CholeskyEliminationTree.h
* @brief   Elimination tree structure to perform symbolic factorization and Cholesky factorization
* @author  Roger Hsiao (rogerhh)
* @date    Feb. 8, 2023
*/

#pragma once

#include <gtsam/base/Matrix.h>
#include <gtsam/base/MatrixSerialization.h>
#include <gtsam/base/FastVector.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/linear/VectorValues.h>
#include <Eigen/Core>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <iostream>
#include <gtsam/base/LowerTriangularColumnMatrix.h>
#include <gtsam/base/SparseLowerTriangularBlockMatrix.h>
#include <gtsam/base/SparseSymmetricBlockMatrix.h>
#include <gtsam/linear/JacobianMatrix.h>
#include <gtsam/base/Workspace.h>

namespace gtsam {

class GTSAM_EXPORT CholeskyEliminationTree {
private:
    class Node;
    class Clique;

    friend std::ostream& operator<<(std::ostream& os, const Clique& clique);

public:
    typedef SparseColumnBlockMatrix::RowMajorMatrix RowMajorMatrix;
    typedef SparseColumnBlockMatrix::ColMajorMatrix ColMajorMatrix;
    typedef SparseColumnBlockMatrix::Block Block;
    typedef SparseColumnBlockMatrix::constBlock constBlock;
    typedef std::pair<size_t, size_t> RowHeightPair;
    typedef LowerTriangularColumnMatrix::BlockIndexVector BlockIndexVector;
    typedef std::shared_ptr<Node> sharedNode;
    typedef std::shared_ptr<Clique> sharedClique;
    typedef std::shared_ptr<NonlinearFactor> sharedFactor;

private:

    std::vector<sharedNode> nodes_;
    sharedClique root_;
    std::vector<size_t> widths_;
    SparseLowerTriangularBlockMatrix cholesky_;
    SparseColumnBlockMatrix delta_;   // \Delta // We're solving AtA Delta = Atb
    std::vector<size_t> ordering_; // Key to index
    std::vector<Key> orderingToKey_; // index to Key
    std::vector<sharedFactor> factors_;
    std::vector<ColMajorMatrix> cachedLinearFactors_;
    Workspace workspace;

    enum MarkedStatus {UNMARKED, EDIT, RECONSTRUCT, NEW};
    std::vector<MarkedStatus> markedStatus_;

    enum LinearizeStatus {UNLINEARIZED, NEWLINEARIZED, RELINEARIZED, LINEARIZED};
    
    std::vector<LinearizeStatus> factorLinearizeStatus_;

    std::vector<std::vector<Key>> descendants_;
    std::vector<std::vector<Key>> changedDescendants_;

    std::vector<bool> backSolveKeys_;
    std::vector<Key> deltaReplaceMask_;
    std::vector<bool> is_reordered_;

    bool groupCliqueFlag = false;
    size_t lowestReorderedIndex = -1;

    static const size_t LAST_ROW;

public:
    CholeskyEliminationTree();

    void addVariables(const Values& newTheta);

    // Mark directly changed keys and keys that we explicitly want to update (extraKeys)
    // observedKeys are the keys associated with the new factors
    void markAffectedKeys(const NonlinearFactorGraph& nonlinearFactors,
                          const FactorIndices& newFactorIndices,
                          const KeySet& relinKeys, 
                          const std::optional<FastList<Key>>& extraKeys,
                          KeySet* affectedKeys,
                          KeySet* observedKeys);

    // Mark all ancestors of directly changed keys and disconnect child from parent 
    // As parent might change. But do not reset colStructure as that does not change
    void markAncestors(const KeySet& affectedKeys, KeySet* markedKeys);

    void symbolicElimination(const KeySet& markedKeys);

    void choleskyElimination(const Values& theta);

    void updateOrdering(const KeySet& markedKeys, const KeySet& observedKeys);
    
    // Go from top of tree and solve for delta. Lt x = y
    // Forward solve is already incorporated in Cholesky
    void backwardSolve(VectorValues* delta_ptr, double tol = 0.001);

    void updateDelta(VectorValues* delta_ptr) const;

    void print(std::ostream& os);

private:
    size_t colWidth(const Key key) const;

    void markKey(const Key key, KeySet* markedKeys);

    void symbolicEliminateKey(const Key key);

    void mergeChildColStructure(sharedNode parent, sharedNode child);

    // Allocate required memory for outer products
    // Required memory is the max sum of matrix sizes for any path from root to leaf
    void allocateStack();
    void deallocateStack();

    struct OrderingLess {
        CholeskyEliminationTree* eTreePtr = nullptr;
        OrderingLess() : eTreePtr(nullptr) {}
        OrderingLess(CholeskyEliminationTree* eTreePtr_in) : eTreePtr(eTreePtr_in) {}
        bool operator()(const Key lhs, const Key rhs) const {
            assert(eTreePtr);
            return eTreePtr->ordering_.at(lhs) < eTreePtr->ordering_.at(rhs);
        }

        bool operator()(const std::pair<Key, RowHeightPair>& lhs, 
                        const std::pair<Key, RowHeightPair>& rhs) const {
            assert(eTreePtr);
            return eTreePtr->ordering_.at(lhs.first) < eTreePtr->ordering_.at(rhs.first);
        }
    };

    OrderingLess orderingLess;
    size_t ordering_version_ = 0;     // When the ordering updates, only nodes that are touched
                                      // need to reorder blocks. nodes that are not touched
                                      // can delay the update

    void constructCSCMatrix(const std::vector<Key>& reorderKeys,
                            const KeySet& observedKeys,
                            int* nEntries,
                            int* nVars,
                            int* nFactors,
                            std::vector<Key>* keyOrdering,
                            std::vector<int>* A,
                            std::vector<int>* p,
                            std::vector<int>* cmember);

    // Get partial ordering of a subset of variables, this subset will be placed
    // At the end of the total ordering regardless of the tree structure
    // cmember is the constraint group.
    void getPartialReordering(const std::vector<Key>& reorderKeys, 
                              const KeySet& observedKeys,
                              std::vector<Key>* partialOrdering);

    // After reordering, all child cliques of reordered nodes
    // need to be detached and reordered, then attached to appropriate parents
    // This needs to be done before symbolic elimination
    void reparentOrphanClique(sharedClique clique);

    // Allocate and gather columns of the clique to a new workspace
    // If reconstruct is true, will not gather but simply reset the columns
    void allocateAndGatherClique(sharedClique clique, bool allocate, bool reconstruct);

    // Copy columns in the workspace back to the cholesky matrix
    void scatterClique(sharedClique clique, bool scatter=true);

    // Reset the first diagWidth columns in the clique's workspace
    void resetCliqueColumns(sharedClique clique);

    // Only handle edit coming from this clique, if clique is edit, also restore 
    // clique columns, but don't subtract AtA block yet
    void editAndRestoreFromClique(sharedClique clique);

    // Returns true if there is a reconstruct column
    // in that case we need to propagate to children
    void eliminateClique(sharedClique clique);

    // Merge outerproduct of a clique into its parent's workspace
    // If parent is nullptr, do nothing
    void mergeWorkspaceClique(sharedClique clique);

    void restoreClique(sharedClique clique);

    void editOrReconstructFromClique(
            sharedClique clique,
            const std::vector<Key>& destCols,
            double sign);

    void setEditOrReconstruct(sharedClique clique);

    void checkEditOrReconstructClique(sharedClique clique, 
                                      MarkedStatus mode, 
                                      std::vector<Key>* destCols);

    // A is k x m, B is k x n, C is m x n
    void simpleAtB(Block A, Block B, Block C, size_t m, size_t k, size_t n, bool add);

    // ADD AtA of factors. Handle edits of AtA here as well
    void prepareEliminateClique(sharedClique clique, const Values& theta);

    // Reorder cliques that were not reset during the reordering
    // Each node needs to have its data structures reordered
    void reorderClique(sharedClique clique);

    void backwardSolveClique(sharedClique clique, VectorValues* delta_ptr, double tol = 0.001);
    
    void checkInvariant_afterMarkAffected() const;
    void checkInvariant_afterMarkAncestors() const;
    void checkInvariant_afterSymbolic() const;
    void checkInvariant_afterCholesky() const;
    void checkInvariant_afterBackSolve() const;

    bool sorted_no_duplicates(const std::vector<Key>& v) const;

};


} // namespace gtsam

