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
#include <gtsam/base/SparseLowerTriangularBlockMatrix.h>
#include <gtsam/base/SparseSymmetricBlockMatrix.h>
#include <gtsam/linear/JacobianMatrix.h>

namespace gtsam {

class GTSAM_EXPORT CholeskyEliminationTree {
private:
    class Node;
    class ColumnStructure;
    typedef SparseColumnBlockMatrix::RowMajorMatrix RowMajorMatrix;
    typedef SparseColumnBlockMatrix::Block Block;
    typedef SparseColumnBlockMatrix::constBlock constBlock;
    typedef std::pair<size_t, size_t> RowHeightPair;
    typedef std::shared_ptr<Node> sharedNode;
    typedef std::shared_ptr<NonlinearFactor> sharedFactor;

    std::vector<sharedNode> nodes_;
    sharedNode root_;
    // SparseSymmetricBlockMatrix hessian_;
    SparseLowerTriangularBlockMatrix cholesky_;
    SparseColumnBlockMatrix b_;   
    SparseColumnBlockMatrix y_;   // L^-1 Atb
    SparseColumnBlockMatrix delta_;   // \Delta // We're solving AtA Delta = Atb
    JacobianMatrix jacobian_;
    std::vector<size_t> ordering_; // Key to index
    std::vector<sharedFactor> factors_;
    // How many nonzero entries in the A matrix
    size_t nEntries_ = 0;


    enum LinearizeStatus {UNLINEARIZED, RELINEARIZED, LINEARIZED};
    
    std::vector<LinearizeStatus> factorLinearizeStatus_;

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

    void choleskyElimination(Values& theta);

    void updateOrdering(KeySet* affectedKeys);
    
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
        bool operator()(const Key lhs, const Key rhs) {
            assert(eTreePtr);
            return eTreePtr->ordering_.at(lhs) < eTreePtr->ordering_.at(rhs);
        }
    };

    OrderingLess orderingLess;

    void mergeColStructure(const Key key,
                           const Key ignoreChildKey,
                           const std::map<Key, size_t, OrderingLess>& src, 
                           std::map<Key, size_t, OrderingLess>* dest);

    void mergeColContribution(const Key key,
                              const Key ignoreChildKey,
                              const std::unordered_map<Key, size_t>& src,
                              std::unordered_map<Key, size_t>* dest);

    // set the new variable ordering. TODO: Add an option for metis
    void getTotalReordering();

    void linearFactor(sharedFactor factor, const Values& theta, 
                      std::vector<RowMajorMatrix>* A,
                      Vector* b);

    // Eliminate the column. The column matrix should already be set up
    void eliminateColumn(SparseColumnBlockMatrix* column_ptr, 
                         Eigen::VectorXd* diagonalY);

    void handleEdits(sharedNode node);

    void handleReconstructs(sharedNode node);

    // Prepares the column for edit. Multiply cholesky of diagonal block to the whole column
    // Then SUBTRACT AtA of factors that need t obe relinearized
    void prepareEditColumn(sharedNode node);


    void backwardSolveNode(sharedNode node);

};

} // namespace gtsam


// #pragma once
// 
// #include <gtsam/base/Matrix.h>
// #include <gtsam/base/MatrixSerialization.h>
// #include <gtsam/base/FastVector.h>
// #include <gtsam/nonlinear/NonlinearFactorGraph.h>
// #include <Eigen/Core>
// #include <vector>
// #include <unordered_map>
// #include <unordered_set>
// #include <utility>
// #include <iostream>
// #include <gtsam/base/SparseUpperTriangularBlockMatrix.h>
// #include <gtsam/base/SparseSymmetricBlockMatrix.h>
// 
// namespace gtsam {
// 
//     class GTSAM_EXPORT CholeskyEliminationTree
//     {
// 
//     private:
//         // Each node represents a variable (TODO: figure out how to make supernodes)
//         struct Node {
//             typedef std::shared_ptr<Node> sharedNode;
//             typedef std::shared_ptr<NonlinearFactor> sharedFactor;
//             Key key;
//             sharedNode parent = nullptr;
//             std::unordered_set<sharedNode> children;
//             std::unordered_set<Key> fillInKeys;         // All keys we interact with that has a higher order than us
//             std::unordered_set<Key> conditionalKeys;    // All keys we interact with that has a lower order than us
//             std::unordered_map<Key, size_t> keyFactorCount;     // Number of times this Key interacts with another Key
//             FactorIndexSet factorIndexSet;
//             bool relinearized = false;
//             bool marked = false;    
// 
//             Node(const Key key_in) 
//             : key(key_in) {}
//             // There are two cases in which we need to recompute the AtA of all factors in a column
//             // 1. A column is relinearized. In this case, all Hessian blocks (i,j) needs to be set to zero where A_ki, A_kj != 0, and the column is set to 0
//             // 2. A factor between (i, j) is removed. In this case, the column block i and j needs to be recomputed. We only need to compute H(i, j), H(i, j), and H(j, j). Actually, we might as well mark the two columns as being relinearized, and bank on this operation not happening a lot. We could also cache the jacobian, but it seems like too much memory requirement for an operation that doesn't happen very often. 
// 
//             // compute the jacobian of a factor given the values v
//             // TODO: Do we need to cache the Jacobian? Should we just cache the Hessian?
//             // The hessian should just be a sum of all AtA matrices from Jacobian blocks
//             // If there is a block deletion or a relinearization
//             void linearizefactor(sharedFactor factor, Values v);
// 
//             // Do AtA for a Jacobian blocks in this factor to variables that have a equal or larger keyIndex than ours
//             // factors to variables that have a smaller keyIndex are processed already
//             // And live in contribution matrices
//             void constructHessian();
// 
//             // do cholesky on self
//             // AtA = RtR
//             // S = R^-1B
//             // C -= StS
//             void CholeskyPartial();
// 
// 
//             // Eliminate this node. All information should be propagated up to parent
//             void eliminate() const {
//                 // linearize all factors in eliminateFactors
//                 // construct hessian in eliminate Factors
//                 // Do Cholesky on self
//                 //
//             }
// 
//             void print(const std::string& str, const KeyFormatter& keyFormatter) const;
//         };
// 
//         typedef std::shared_ptr<Node> sharedNode;
//         typedef std::shared_ptr<NonlinearFactor> sharedFactor;
// 
//     public:
//         void addVariables(const Values& newTheta); 
// 
//         void updateFactorsAndMarkAffectedKeys(const NonlinearFactorGraph& nonlinearFactors,
//                                               const FactorIndices& newFactorIndices, 
//                                               const FactorIndices& removeFactorIndices,
//                                               const std::optional<FastList<Key>>& extraKeys,
//                                               KeySet& affectedKeys);
// 
// 
//         // For any observed keys and relin key, traverse up the tree and find all their ancestors
//         // Mark all such nodes, and detach them from their parents
//         // Then take all the children of the marked keys, and detach them
//         void markAllAffectedKeys(const KeySet& observedKeys, 
//                                  const KeySet& relinKeys,
//                                  KeySet* markedKeys,
//                                  KeyVector* orphanKeys);
// 
//         // Given some ordering, reintegrate all the nodes back into the tree
//         // While doing that, collect the size each column of the R matrix needs to allocate 
//         // and allocate it at the end
//         // We can resize all the column blocks at relinearization
//         void symbolicElimination(const Ordering& ordering, const KeyVector& orphanKeys);
// 
//         void choleskyElimination();
// 
//         void printTreeStructure(std::ostream& os);
// 
// 
//     private:
//         // Mark this node for re-elimination
//         void markNode(sharedNode node, KeySet* markedKeys);
// 
//         // Add all children of this node to a set of orphan keys
//         // And disconnect them from the parent
//         void orphanChildren(sharedNode node, KeyVector* orphanKeys);
// 
//         // 0. Assume all contribution blocks are set
//         // 1. Assume all factors are relinearized. Note: There should be a more cache friendly way to relinearize factors while keeping track of old contribution blocks
//         // 2. Take dense Cholesky of diagonal block
//         // 3. Solve all blocks under the diagonal
//         // 4. Compute fill-in blocks
//         void choleskyEliminateNode(sharedNode node);
// 
//         void reParentOrphan(sharedNode node, const std::unordered_map<Key, size_t>& keyToOrdering);
//         void resolveAllocate();
// 
//         void addFillInKey(sharedNode node, const Key otherKey);
// 
//         void symbolicEliminateNode(sharedNode node,
//                                    size_t myOrder,
//                                    const std::unordered_map<Key, size_t>& keyToOrdering);
//         
//         // Mark a single relinKey
//         // And then we need to mark the ancestors of those keys
//         void markRelinKey(const Key relinKey);
// 
//         void validateTree();
// 
//         std::vector<sharedNode> nodes_; 
//         sharedNode root_;
//         std::unordered_map<FactorIndex, bool> factorLinearizeStatus_;
// 
//         SparseUpperTriangularBlockMatrix cholesky_;
//         SparseSymmetricBlockMatrix hessian_;
// 
//         // When reparenting, we need to add orphan fillins to parent
//         // Actually need to reparent first
// 
//     };
// }
