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

from copy import deepcopy

from utils.utils import *

class ConvergenceChecker:
    # atol is the absolute graph error, rtol is how much the graph error reduced
    def __init__(self, max_iter=10, atol=0.0017, rtol=0.01):
        self.max_iter = max_iter
        self.iter = 0
        self.atol = atol
        self.rtol = rtol
        self.last_chi2 = None

    def check_convergence(self, graph, estimate, factor_dim):

        chi2_graph = chi2_red(graph, estimate, factor_dim)

        if chi2_graph < self.atol:
            return True, chi2_graph

        if self.last_chi2 is not None:
            diff = math.fabs(chi2_graph - self.last_chi2)
            print("r error = ", diff / chi2_graph)
            if diff / chi2_graph < self.rtol:
                return True, chi2_graph

        self.last_chi2 = chi2_graph

        self.iter += 1
        if self.iter >= self.max_iter:
            return True, chi2_graph

        return False, chi2_graph
