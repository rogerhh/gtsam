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

class DirectSolver:
    def __init__(self):
        pass

    # Solve an Mx = y problem directly. M has to be SPD
    def solve(self, M, y):
        chol_factor = cholesky(M)
        x = chol_factor(y)
        return x.toarray()
