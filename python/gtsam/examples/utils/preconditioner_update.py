"""
preconditioner_update.py: Given an old preconditioner, an old A and a new A, return a new preconditioner
"""

from sksparse.cholmod import cholesky, cholesky_AAt
from scipy.sparse import csc_matrix, csr_matrix
from scipy.sparse.linalg import cg, spsolve_triangular, inv
import matplotlib.pyplot as plt
import scipy
import numpy as np
from copy import deepcopy

from utils.utils import permute, unpermute

import gtsam
from functools import cmp_to_key

class TrueIdentityPreconditionerUpdater:
    def updatePreconditioner(self, L_old=None, P_old=None, P_new=None, \
                                   Lamb_old=None, Lamb_new=None, \
                                   A_old=None, A_prime=None, A_tilde=None):
        height_old, width_old = Lamb_old.shape
        height_new, width_new = Lamb_new.shape
        return scipy.sparse.eye(height_new)

class IdentityPreconditionerUpdater:
    def updatePreconditioner(self, L_old=None, P_old=None, P_new=None, \
                                   Lamb_old=None, Lamb_new=None, \
                                   A_old=None, A_prime=None, A_tilde=None):
        height_old, width_old = Lamb_old.shape
        height_new, width_new = Lamb_new.shape
        L_new = deepcopy(L_old)
        L_new.resize(Lamb_new.shape)

        diff = height_new - height_old
        new_diag = L_new.diagonal()
        new_diag[height_old:height_new] = np.ones((diff,))
        L_new.setdiag(new_diag)

        return L_new

class ExtendedDiagonalPreconditionerUpdater:
    def updatePreconditioner(self, L_old=None, P_old=None, P_new=None, \
                                   Lamb_old=None, Lamb_new=None, \
                                   A_old=None, A_prime=None, A_tilde=None):
        height_old, width_old = Lamb_old.shape
        height_new, width_new = Lamb_new.shape
        L_new = deepcopy(L_old)
        L_new.resize(Lamb_new.shape)

        diff = height_new - height_old
        new_diag = L_new.diagonal()
        new_diag[height_old:height_new] = np.sqrt(Lamb_new[height_old:height_new, height_old:height_new].diagonal())
        L_new.setdiag(new_diag)

        return L_new

class IncompleteCholeskyOnDiagPreconditionerUpdater:
    def updatePreconditioner(self, L_old=None, P_old=None, P_new=None, \
                                   Lamb_old=None, Lamb_new=None, \
                                   A_old=None, A_prime=None, A_tilde=None):
        height_old, width_old = Lamb_old.shape
        height_new, width_new = Lamb_new.shape
        L_new = deepcopy(L_old)
        L_new.resize(Lamb_new.shape)

        diff = height_new - height_old
        new_diag = L_new.diagonal()
        new_diag[height_old:height_new] = np.sqrt(Lamb_new[height_old:height_new, height_old:height_new].diagonal())
        L_new.setdiag(new_diag)

        return L_new

class CholeskyOnNewBlockPreconditionerUpdater:
    def updatePreconditioner(self, L_old=None, P_old=None, P_new=None, \
                                   Lamb_old=None, Lamb_new=None, \
                                   A_old=None, A_prime=None, A_tilde=None):
        height_old, width_old = Lamb_old.shape
        height_new, width_new = Lamb_new.shape
        L_new = deepcopy(L_old)
        L_new.resize(Lamb_new.shape)

        new_block = Lamb_new[height_old:height_new,height_old:height_new]
        chol_block = cholesky(new_block, ordering_method="natural").L()

        L_new[height_old:height_new,height_old:height_new] = chol_block

        return L_new

class IncompleteCholeskyPreconditionerUpdater:
    def updatePreconditioner(self, L_old=None, P_old=None, P_new=None, \
                                   Lamb_old=None, Lamb_new=None, \
                                   A_old=None, A_prime=None, A_tilde=None):
        height_old, width_old = Lamb_old.shape
        height_new, width_new = Lamb_new.shape
        L_new = deepcopy(L_old)
        L_new.resize(Lamb_new.shape)

        Lamb_prime = deepcopy(Lamb_old)
        Lamb_prime.resize(Lamb_new.shape)
        Lamb_prime += A_prime.T @ A_prime

        Lamb_prime_permuted = permute(Lamb_prime, P_new)
        Lamb_prime_old_block = Lamb_prime_permuted[:height_old,:width_old]
        L_new_old_block = cholesky(Lamb_prime_old_block, ordering_method="natural").L()

        # Do the chol update to the complete cholesky
        L_new[L_new.nonzero()] *= 0
        L_new[L_new_old_block.nonzero()] = L_new_old_block[L_new_old_block.nonzero()]

        assert(0, "This code is wrong. See FullIncompleteCholeskyPreconditionerUpdater")
        # Do the chol update to the incomplete cholesky
        for r1, c in zip(*Lamb_prime_permuted[height_old:height_new, :].nonzero()):
            r = r1 + height_old
            if c > r:
                continue
            assert(L_new[r, c] == 0)

            diag_val = Lamb_prime_permuted[c, c]
            contri_val = (L_new[0:c-1, r].T @ L_new[0:c-1, c])[0, 0]
            if c == r:
                L_new[c, c] = np.sqrt(diag_val - contri_val)
            else:
                subdiag_val = Lamb_prime_permuted[r, c]
                L_new[r, c] = (subdiag_val - contri_val) / diag_val

        return L_new

class IncompleteCholeskyWithRelinLambdaPreconditionerUpdater:
    def updatePreconditioner(self, L_old=None, P_old=None, P_new=None, \
                                   Lamb_old=None, Lamb_new=None, \
                                   A_old=None, A_prime=None, A_tilde=None):
        height_old, width_old = Lamb_old.shape
        height_new, width_new = Lamb_new.shape
        L_new = deepcopy(L_old)
        L_new.resize(Lamb_new.shape)

        Lamb_new_permuted = permute(Lamb_new, P_new)
        Lamb_prime_old_block = Lamb_new_permuted[:height_old,:width_old]
        L_new_old_block = cholesky(Lamb_prime_old_block, ordering_method="natural").L()

        # Do the chol update to the complete cholesky
        L_new[L_new.nonzero()] *= 0
        L_new[L_new_old_block.nonzero()] = L_new_old_block[L_new_old_block.nonzero()]

        # Do the chol update to the incomplete cholesky
        for r1, c in zip(*Lamb_new_permuted[height_old:height_new, :].nonzero()):
            r = r1 + height_old
            if c > r:
                continue
            assert(L_new[r, c] == 0)

            diag_val = Lamb_new_permuted[c, c]
            contri_val = (L_new[0:c-1, r].T @ L_new[0:c-1, c])[0, 0]
            if c == r:
                L_new[c, c] = np.sqrt(diag_val - contri_val)
            else:
                subdiag_val = Lamb_new_permuted[r, c]
                L_new[r, c] = (subdiag_val - contri_val) / diag_val

        assert(0, "This code is wrong. See FullIncompleteCholeskyPreconditionerUpdater")
        return L_new

class FullIncompleteCholeskyPreconditionerUpdater:
    # No more complete Cholesky. Only incomplete cholesky for all high norm rows
    def updatePreconditioner(self, L_old=None, P_old=None, P_new=None, \
                                   Lamb_old=None, Lamb_new=None, \
                                   A_old=None, A_prime=None, A_tilde=None):
        height_old, width_old = Lamb_old.shape
        height_new, width_new = Lamb_new.shape
        L_new = deepcopy(L_old)
        L_new.resize(Lamb_new.shape)

        Lamb_prime = deepcopy(Lamb_old)
        Lamb_prime.resize(Lamb_new.shape)
        Lamb_prime += A_prime.T @ A_prime

        # Lamb_prime_permuted = permute(Lamb_prime, P_new)
        Lamb_prime_permuted = permute(Lamb_new, P_new)

        L_ref = cholesky(Lamb_prime_permuted, ordering_method='natural').L()

        L_new = deepcopy(Lamb_prime_permuted)
        L_new[L_new.nonzero()] *= 0

        rows, cols = Lamb_prime_permuted.nonzero()

        def comp(i, j):
            if cols[i] < cols[j]:
                return -1
            elif cols[i] > cols[j]:
                return 1
            else:
                return rows[i] - rows[j]

        indices = sorted(range(len(rows)), key=cmp_to_key(comp))
        cols = cols[indices]
        rows = rows[indices]

        # Do the chol update to the incomplete cholesky
        for r, c in zip(rows, cols):
            if c > r:
                continue
            assert(L_new[r, c] == 0)

            contri_val = 0
            if c > 0:
                contri_val = (L_new[c, 0:c] @ L_new[r, 0:c].T)[0, 0]
            if c == r:
                L_new[c, c] = np.sqrt(Lamb_prime_permuted[c, c] - contri_val)
            else:
                L_new[r, c] = (Lamb_prime_permuted[r, c] - contri_val) / L_new[c, c]

        L_new.eliminate_zeros()

        return L_new


class CholeskyUpdatePreconditionerUpdater:
    def updatePreconditioner(self, L_old=None, P_old=None, P_new=None, \
                                   Lamb_old=None, Lamb_new=None, \
                                   A_old=None, A_prime=None, A_tilde=None):
        height_old, width_old = Lamb_old.shape
        height_new, width_new = Lamb_new.shape
        L_new = deepcopy(L_old)
        L_new.resize(Lamb_new.shape)

        Lamb_prime = deepcopy(Lamb_old)
        Lamb_prime.resize(Lamb_new.shape)
        Lamb_prime += A_prime.T @ A_prime

        Lamb_prime_permuted = permute(Lamb_prime, P_new)
        L_new = cholesky(Lamb_prime_permuted, ordering_method="natural").L()

        return L_new

# This is where we take Cholesky update of the high norm rows
class SelectiveCholeskyUpdatePreconditionerUpdater:
    def updatePreconditioner(self, L_old=None, P_old=None, P_new=None, \
                                   Lamb_old=None, Lamb_new=None, \
                                   A_old=None, A_prime=None, A_tilde=None):

        A_prime_height, A_prime_width = A_prime.shape
        A_tilde_height, A_tilde_width = A_tilde.shape

        sel_thresh = 0.1
        factor_dim = 6
        high_rows = np.where(abs(A_tilde.A) > sel_thresh)[0]
        high_rows = np.unique(high_rows)

        print("length of A_tilde_high = ", A_tilde[high_rows].shape)
        A_tilde_low = deepcopy(A_tilde)
        A_tilde_low[high_rows] = 0
        A_tilde_high = A_tilde - A_tilde_low

        height_old, width_old = Lamb_old.shape
        height_new, width_new = Lamb_new.shape

        A_prime_tilde_high = deepcopy(A_old)
        A_prime_tilde_high.resize((A_tilde_height + A_prime_height, A_prime_width))

        A_prime_tilde_high[:A_tilde_height, :A_tilde_width] += A_tilde_high
        A_prime_tilde_high[A_tilde_height:, :] += A_prime
        Lamb_prime = A_prime_tilde_high.T @ A_prime_tilde_high

        Lamb_prime_permuted = permute(Lamb_prime, P_new)
        L_new = cholesky(Lamb_prime_permuted, ordering_method="natural").L()

        return L_new

# This is where we take Cholesky update of the high norm row-blocks 
class SelectiveCholeskyUpdatePreconditionerUpdater2:
    def updatePreconditioner(self, L_old=None, P_old=None, P_new=None, \
                                   Lamb_old=None, Lamb_new=None, \
                                   A_old=None, A_prime=None, A_tilde=None):

        A_prime_height, A_prime_width = A_prime.shape
        A_tilde_height, A_tilde_width = A_tilde.shape

        sel_thresh = 0.1
        factor_dim = 6
        high_factors = np.where(abs(A_tilde.A) > sel_thresh)[0] // factor_dim
        high_factors = np.unique(high_factors)
        high_rows = []
        for factor in high_factors:
            high_rows.extend(range(factor * factor_dim, (factor + 1) * factor_dim))

        print("length of A_tilde_high = ", A_tilde[high_rows].shape)
        A_tilde_low = deepcopy(A_tilde)
        A_tilde_low[high_rows] = 0
        A_tilde_high = A_tilde - A_tilde_low

        height_old, width_old = Lamb_old.shape
        height_new, width_new = Lamb_new.shape

        A_prime_tilde_high = deepcopy(A_old)
        A_prime_tilde_high.resize((A_tilde_height + A_prime_height, A_prime_width))

        A_prime_tilde_high[:A_tilde_height, :A_tilde_width] += A_tilde_high
        A_prime_tilde_high[A_tilde_height:, :] += A_prime
        Lamb_prime = A_prime_tilde_high.T @ A_prime_tilde_high

        Lamb_prime_permuted = permute(Lamb_prime, P_new)
        L_new = cholesky(Lamb_prime_permuted, ordering_method="natural").L()

        return L_new

# Incomplete Cholesky with nnz(col) + k nonzeros per column
# The k new nonzeros are the k highest nonzeros in the column
# Update incomplete Cholesky factorization with high norm rows
class IncompleteCholeskyPreconditionerUpdater2:
    def __init__(self, k):
        self.k = k

    def updatePreconditioner(self, L_old=None, P_old=None, P_new=None, \
                                   Lamb_old=None, Lamb_new=None, \
                                   A_old=None, A_prime=None, A_tilde=None):

        A_prime_height, A_prime_width = A_prime.shape
        A_tilde_height, A_tilde_width = A_tilde.shape

        sel_thresh = 0.1
        factor_dim = 6
        high_factors = np.where(abs(A_tilde.A) > sel_thresh)[0] // factor_dim
        high_factors = np.unique(high_factors)
        high_rows = []
        for factor in high_factors:
            high_rows.extend(range(factor * factor_dim, (factor + 1) * factor_dim))

        print("length of A_tilde_high = ", A_tilde[high_rows].shape)
        A_tilde_low = deepcopy(A_tilde)
        A_tilde_low[high_rows] = 0
        A_tilde_high = A_tilde - A_tilde_low

        height_old, width_old = Lamb_old.shape
        height_new, width_new = Lamb_new.shape

        A_prime_tilde_high = deepcopy(A_old)
        A_prime_tilde_high.resize((A_tilde_height + A_prime_height, A_prime_width))

        A_prime_tilde_high[:A_tilde_height, :A_tilde_width] += A_tilde_high
        A_prime_tilde_high[A_tilde_height:, :] += A_prime
        Lamb_prime = A_prime_tilde_high.T @ A_prime_tilde_high

        Lamb_prime_permuted = permute(Lamb_prime, P_new)
        L_new = cholesky(Lamb_prime_permuted, ordering_method="natural").L()

        return L_new
        height_old, width_old = Lamb_old.shape
        height_new, width_new = Lamb_new.shape
        L_new = deepcopy(L_old)
        L_new.resize(Lamb_new.shape)

        Lamb_prime = deepcopy(Lamb_old)
        Lamb_prime.resize(Lamb_new.shape)
        Lamb_prime += A_prime.T @ A_prime

        # Lamb_prime_permuted = permute(Lamb_prime, P_new)
        Lamb_prime_permuted = permute(Lamb_new, P_new)

        L_ref = cholesky(Lamb_prime_permuted, ordering_method='natural').L()

        L_new = deepcopy(Lamb_prime_permuted)
        L_new[L_new.nonzero()] *= 0

        rows, cols = Lamb_prime_permuted.nonzero()

        def comp(i, j):
            if cols[i] < cols[j]:
                return -1
            elif cols[i] > cols[j]:
                return 1
            else:
                return rows[i] - rows[j]

        indices = sorted(range(len(rows)), key=cmp_to_key(comp))
        cols = cols[indices]
        rows = rows[indices]

        # Do the chol update to the incomplete cholesky
        for r, c in zip(rows, cols):
            if c > r:
                continue
            assert(L_new[r, c] == 0)

            contri_val = 0
            if c > 0:
                contri_val = (L_new[c, 0:c] @ L_new[r, 0:c].T)[0, 0]
            if c == r:
                L_new[c, c] = np.sqrt(Lamb_prime_permuted[c, c] - contri_val)
            else:
                L_new[r, c] = (Lamb_prime_permuted[r, c] - contri_val) / L_new[c, c]

        L_new.eliminate_zeros()

        return L_new
