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

namespace gtsam {

class GTSAM_EXPORT CholeskyEliminationTree {
private:
    class Node;
    class Clique;

public:
    typedef SparseColumnBlockMatrix::RowMajorMatrix RowMajorMatrix;
    typedef SparseColumnBlockMatrix::ColMajorMatrix ColMajorMatrix;
    typedef SparseColumnBlockMatrix::Block Block;
    typedef SparseColumnBlockMatrix::constBlock constBlock;
    typedef std::pair<size_t, size_t> RowHeightPair;
    typedef std::shared_ptr<Node> sharedNode;
    typedef std::shared_ptr<Clique> sharedClique;
    typedef std::shared_ptr<NonlinearFactor> sharedFactor;

private:

    std::vector<sharedNode> nodes_;
    // sharedNode root_;
    sharedClique root_;
    // SparseSymmetricBlockMatrix hessian_;
    SparseLowerTriangularBlockMatrix cholesky_;
    // SparseColumnBlockMatrix b_;  // This is incorporated in the Jacobian
    // SparseColumnBlockMatrix y_;   // L^-1 Atb // This is incorporated in choleksy
    SparseColumnBlockMatrix delta_;   // \Delta // We're solving AtA Delta = Atb
    JacobianMatrix jacobian_;
    std::vector<size_t> ordering_; // Key to index
    std::vector<Key> orderingToKey; // index to Key
    std::vector<sharedFactor> factors_;
    // How many nonzero entries in the A matrix
    size_t nEntries_ = 0;

    enum LinearizeStatus {UNLINEARIZED, NEWLINEARIZED, RELINEARIZED, LINEARIZED};
    
    std::vector<LinearizeStatus> factorLinearizeStatus_;

    static const size_t LAST_ROW = LowerTriangularColumnMatrix::LAST_ROW;

public:
    CholeskyEliminationTree();

    void addVariables(const Values& newTheta);

    // Mark directly changed keys and keys that we explicitly want to update (extraKeys)
    void markAffectedKeys(const NonlinearFactorGraph& nonlinearFactors,
                          const FactorIndices& newFactorIndices,
                          const KeySet& relinKeys, 
                          const std::optional<FastList<Key>>& extraKeys,
                          KeySet* affectedKeys);

    // Mark all ancestors of directly changed keys and disconnect child from parent 
    // As parent might change. But do not reset colStructure as that does not change
    void markAncestors(const KeySet& affectedKeys, KeySet* markedKeys);

    // Mark all ancestors of directly changed keys, 
    void markAncestorsForReordering(const KeySet& affectedKeys, KeySet* markedKeys);

    // There should be two versions of symbolic elim
    // One with reordering and one without
    void symbolicEliminationWithOrdering(const Ordering& ordering);

    void symbolicElimination(const KeySet& markedKeys);

    void choleskyElimination(const Values& theta);


    // void updateOrdering(KeySet* affectedKeys);
    void updateOrdering(const KeySet& markedKeys);
    
    // Go from top of tree and solve for delta. Lt x = y
    // Forward solve is already incorporated in Cholesky
    void backwardSolve();

    void updateDelta(VectorValues* delta_ptr) const;

    void print(std::ostream& os);

private:
    void markKey(const Key key, KeySet* markedKeys);

    void symbolicEliminateKey(const Key key);

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

    void mergeColStructure(const Key key,
                           const Key ignoreChildKey,
                           const std::map<Key, size_t, OrderingLess>& src, 
                           std::map<Key, size_t, OrderingLess>* dest);

    void mergeColContribution(const Key key,
                              const Key ignoreChildKey,
                              const std::unordered_map<Key, size_t>& src,
                              std::unordered_map<Key, size_t>* dest);

    void constructCSCMatrix(const std::vector<Key>& reorderKeys,
                            int* nEntries,
                            int* nVars,
                            int* nFactors,
                            std::vector<int>* A,
                            std::vector<int>* p,
                            std::vector<int>* cmember,
                            std::unordered_set<Key>* is_reordered);

    // set the new variable ordering. TODO: Add an option for metis
    void getTotalReordering();

    // Get partial ordering of a subset of variables, this subset will be placed
    // At the end of the total ordering regardless of the tree structure
    // cmember is the constraint group.
    void getPartialReordering(const std::vector<Key>& reorderKeys, 
                              std::vector<Key>* partialOrdering);

    // After reordering, all child cliques of reordered nodes
    // need to be detached and reordered, then attached to appropriate parents
    // This needs to be done before symbolic elimination
    void reparentOrphanClique(sharedClique clique);

    std::vector<sharedClique> orphanCliques;

    void linearFactor(sharedFactor factor, const Values& theta, 
                      std::vector<RowMajorMatrix>* A,
                      Vector* b);

    // Eliminate the column. The column matrix should already be set up
    // void eliminateColumn(SparseColumnBlockMatrix* column_ptr, 
    //                      Eigen::VectorXd* diagonalY);

    // Only gather columns that are allocated
    // If nothing to gather (e.g., column might be new) then return false
    bool gatherColumns(sharedClique clique, 
                       ColMajorMatrix* m,
                       size_t* totalWidth,
                       size_t* totalHeight,
                       std::vector<Key>* keys,
                       std::vector<size_t>* widths,
                       std::vector<size_t>* heights,
                       std::vector<std::pair<Key, RowHeightPair>>* blockStartVec,
                       std::vector<LowerTriangularColumnMatrix*>* columns);

    void scatterColumns(const ColMajorMatrix& m, 
                        const std::vector<size_t>& widths,
                        const std::vector<size_t>& heights,
                        const std::vector<LowerTriangularColumnMatrix*>& columns);

    // Only handle edit coming from this clique, if clique is edit, also restore 
    // clique columns, but don't subtract AtA block yet
    void editFromClique(sharedClique clique);
    void reconstructFromClique(sharedClique clique);
    void eliminateClique(sharedClique clique);

    // Handle all edits from this clumn
    // Additional, if this clique is reconstruct, restore column to linear state
    void editAndRestoreFromClique(sharedClique clique);

    void restoreColumn(ColMajorMatrix& m, size_t totalWidth, size_t totalHeight);

    // Given some source columns (children columns) handle edits for keys in editCols
    // Edit keys may not be in the diagonal block
    // m will be changed for the edited columns
    // editCols should only contain keys that the column owns
    void editFromColumn(const ColMajorMatrix& m, 
            const std::vector<std::pair<Key, RowHeightPair>>& blockStartVec,
            const std::vector<Key>& editCols);
    void handleEdits(sharedNode node);

    // Given some source columns, handle reconstructs for keys in reconstructCols
    // m will not be changed as this supercolumn must be unmarked
    void reconstructFromColumn(const ColMajorMatrix& m, 
            const std::vector<std::pair<Key, RowHeightPair>>& blockStartVec,
            const std::vector<Key>& reconstructKeys);

    void eliminateColumn(ColMajorMatrix& m,
            const size_t totalWidth,
            const size_t totalHeight,
            const std::vector<std::pair<Key, RowHeightPair>>& blockStartVec);

    // Return a destBlockStart that constains all the rows and heights of the destination
    // blocks in the destination column, merging blocks when contiguous
    template<typename Iterator>
    void findDestBlocks(Iterator start, Iterator end, 
                        const std::unordered_map<Key, RowHeightPair>& destBlockStartMap,
                        std::vector<RowHeightPair>* destBlockStart);

    void setEditOrReconstruct(sharedClique clique);

    void handleReconstructs(sharedNode node);

    // Prepares the column for edit.
    // SUBTRACT AtA of factors that need to be relinearized
    void prepareEditClique(sharedClique clique);
    void prepareEditColumn(sharedNode node);

    // ADD AtA of factors
    void prepareEliminateClique(sharedClique clique, const Values& theta);
    void prepareEliminateColumn(sharedNode node, bool is_reconstruct, const Values& theta);

    // Reorder cliques that were not reset during the reordering
    // Each node needs to have its data structures reordered
    void reorderClique(sharedClique clique);

    void prepareReconstructColumn(sharedNode node);

    void backwardSolveColumn(ColMajorMatrix& m,
            const std::vector<std::pair<Key, RowHeightPair>>& blockStartVec);

    void backwardSolveClique(sharedClique clique);
    void backwardSolveNode(sharedNode node);
    
    // Merge clique2 into clique1
    void mergeCliques(sharedClique clique1, sharedClique clique2);

};

} // namespace gtsam

