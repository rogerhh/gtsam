"""
preconditioner_update.py: Given an old preconditioner, an old A and a new A, return a new preconditioner
"""

from scikits.sparse.cholmod import cholesky, cholesky_AAt
from scipy.sparse import csc_matrix, csr_matrix
from scipy.sparse.linalg import cg, spsolve_triangular, inv
import matplotlib.pyplot as plt
import scipy
import numpy as np
from copy import deepcopy

import gtsam

class IdentityPreconditionerUpdater:
    def updatePreconditioner(self, L_old=None, P_old=None, Lamb_old=None, Lamb_new=None, \
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
    def updatePreconditioner(self, L_old=None, P_old=None, Lamb_old=None, Lamb_new=None, \
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
    def updatePreconditioner(self, L_old=None, P_old=None, Lamb_old=None, Lamb_new=None, \
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
