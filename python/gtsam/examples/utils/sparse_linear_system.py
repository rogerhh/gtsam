"""
sparse_linear_system.py: We have a collection of Factor's, where each Factor owns a set of continuous indices. Linearizing/Relinearizing a Factor involves changing the entries of A_data, b_data at the indices that it owns A_data[indices], b_data[indices]. The row and column indices are denoted by A_rows[indices] and A_cols[indices]
A sparse_linear_system owns a list of Factor's, (A_data, A_rows, A_cols) which are 3 vectors that represent a sparse matrix A, and a dense vector b. b can be represented by only 1 vector because all indices are filled and there are no duplicated indices
"""

from scikits.sparse.cholmod import cholesky, cholesky_AAt
from scipy.sparse import csc_matrix, csr_matrix
from scipy.sparse.linalg import cg, spsolve_triangular, inv
import matplotlib.pyplot as plt
import scipy
import numpy as np

import gtsam

class Factor:
    def __init__(self, factor, system):
        self.indices = []
        # A reference to the system so we can change the matrix entries
        self.system = system
        pass

    # Update the entries at the indices that we own
    def linearize(self, theta):
        pass

class SparseLinearSystem:

    def __init__(self, measurements):
        self.factors = []
        self.A_start_rows = []
        self.A_start_indices = []
        self.b_start_indices = []
        self.A_data = []
        self.A_rows = []
        self.A_cols = []
        self.b_data = []
        self.max_rows = 0
        self.max_cols = 0
        self.key_to_col = []
        self.measurements = measurements
        self.measurement_index = 0

    def addVariables(self, new_theta):
        # Only need to add dimension of new theta here
        # Currently only works for Pose2
        for key in new_theta.keys():
            value = new_theta.atPose2(key)
            assert(key == len(self.key_to_col))
            self.key_to_col.append((self.max_cols, 3))
            self.max_cols += 3

    def addFactor(self, factor):
        self.factors.append(factor)
        self.A_start_indices.append(len(self.A_data))
        self.b_start_indices.append(len(self.b_data))

        for key in factor.keys():
            (col, width) = self.key_to_col[key]
            for j in range(width):
                c = col + j
                for i in range(factor.dim()):
                    r = self.max_rows + i
                    self.A_data.append(0)
                    self.A_cols.append(c)
                    self.A_rows.append(r)

        self.max_rows += factor.dim()

        for i in range(factor.dim()):
            self.b_data.append(0)

    # We already own all the factors, just need to add to list of active factors
    def addFactors(self, max_measurement_index):
        for measurement_index in range(self.measurement_index, max_measurement_index):
            self.addFactor(self.measurements.at(measurement_index))
        self.measurement_index = max_measurement_index

    def linearizeFactor(self, factor_index, theta):
        factor = self.factors[factor_index]
        jacobian_factor = factor.linearize(theta)
        matrixA = jacobian_factor.getA()
        vectorb = jacobian_factor.getb()

        A_index = self.A_start_indices[factor_index]
        A_start_row = self.b_start_indices[factor_index]
        col_start = 0

        for key in factor.keys():
            (col, width) = self.key_to_col[key]

            for j in range(width):
                c = col + j
                for i in range(factor.dim()):
                    r = A_start_row + i

                    assert(self.A_rows[A_index] == r)
                    assert(self.A_cols[A_index] == c)

                    self.A_data[A_index] = matrixA[i, col_start + j]

                    A_index += 1

            col_start += width

        b_index = self.b_start_indices[factor_index]
        for i in range(factor.dim()):
            self.b_data[b_index + i] = vectorb[i]

    def linearizeAll(self, theta):
        for factor_index in range(len(self.factors)):
            self.linearizeFactor(factor_index, theta)

    def linearizeNew(self, index_start, index_end, theta):
        for factor_index in range(index_start, index_end):
            self.linearizeFactor(factor_index, theta)

    def linearizeSet(self, relin_factors, theta):
        for factor_index in relin_factors:
            self.linearizeFactor(factor_index, theta)


    # Return A matrix and b vector obtained at step and with linearization point theta
    # Adding all the measurements until step
    def getSystem(self):
        A = csr_matrix((self.A_data, (self.A_rows, self.A_cols)))
        b_rows = range(len(self.b_data))
        b_cols = [0 for _ in b_rows]
        b = csr_matrix((self.b_data, (b_rows, b_cols)))
        # b = np.array(self.b_data).reshape((-1, 1))
        return A, b
