import math

import matplotlib.pyplot as plt
import numpy as np

import gtsam
import gtsam.utils.plot as gtsam_plot

from typing import List

from optparse import OptionParser

from scikits.sparse.cholmod import cholesky, cholesky_AAt
from scipy.sparse import csc_matrix, csr_matrix
from scipy.sparse.linalg import cg, spsolve_triangular, inv
import matplotlib.pyplot as plt
import scipy

from utils.preconditioner_generator import PreconditionerGenerator
from utils.sparse_linear_system import SparseLinearSystem
from utils.problem_advancer import ProblemAdvancer

from copy import deepcopy

cg_count = 0
relin_threshold = 0.001
def cg_count_increment(xk):
    global cg_count
    cg_count += 1

    # relin_count = 0
    # for i in range(0, len(xk), 3):
    #     if np.linalg.norm(xk[i:i+3], ord=np.inf) > relin_threshold:
    #         relin_count += 1
        
    # print(cg_count, relin_count)

class IterativeSolver:
    # L is the approximate cholesky factor used as a preconditioner
    # P is the permutation vector
    def __init__(self, L, P, tol=1e-10):
        self.L = deepcopy(L)
        self.P = deepcopy(P)
        self.cg_tol = tol
        pass

    def apply_preconditioner(self, M, R):
        # RM = spsolve_triangular(R, M, lower=True)
        # RTRM = spsolve_triangular(R.T, RM, lower=False)
        # return RTRM
        return inv(R @ R.T) @ M

    # Solve an Mx = y problem iteratively. M has to be SPD
    # M is assumed to be unpermuted
    def solve(self, M, y, x0=None, old_M=None):

        if x0 is None:
            x0 = np.zeros_like(y)

        # First extend L and P to fit the new matrix
        old_width, _ = self.L.shape
        new_width, _ = M.shape

        L = deepcopy(self.L)
        L.resize((new_width, new_width))

        diag = M.diagonal()
        diag_root = np.sqrt(diag)[old_width:new_width]
        indices = range(old_width, new_width)
        diag_csc = csc_matrix((diag_root, (indices, indices)), shape=(new_width, new_width))
        L += diag_csc

        self.P.resize((new_width,))
        self.P[old_width:new_width] = np.array(range(old_width, new_width))

        # Permute M and y
        P = self.P
        PM = M.toarray()[P[:, np.newaxis], P[np.newaxis, :]]
        Py = y[P]

        # conditioned_M = self.apply_preconditioner(PM, L)
        # conditioned_y = self.apply_preconditioner(Py, L)

        # assert(np.allclose(L @ L.T @ conditioned_M, PM))
        # assert(np.allclose(L @ L.T @ conditioned_y, Py))

        # print(conditioned_M)
        # print(conditioned_y)

        x0 = x0[P]

        global cg_count
        cg_count = 0
        # print("condition number before conditioning: ", np.linalg.cond(PM))
        # print("condition number after conditioning: ", np.linalg.cond((inv(L @ L.T) @ PM)))
        Px, info = cg(PM, Py, x0=x0, M=inv(L @ L.T), callback=cg_count_increment, tol=self.cg_tol)
        # Px, info = cg(PM, Py, x0=x0, callback=cg_count_increment, tol=self.cg_tol)

        print("cg iterations: ", cg_count, info)

        Px = Px.reshape((-1, 1))

        x = np.zeros_like(Px)
        x[P] = Px

        # x, info = cg(M, y, callback=cg_count_increment, tol=self.cg_tol)
        # x = x.reshape((-1, 1))
        # print("cg iterations: ", cg_count, info)

        return x
