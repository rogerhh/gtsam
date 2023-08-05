"""
linear_operator.py: A subclass of scipy.linalg.LinearOperator. Given A, L, implements matvec(v) that returns L^-1A^TAL^-Tv, but applies each matrix to v sequentially
"""

from scipy.sparse.linalg import LinearOperator
from scipy.sparse.linalg import cg, spsolve_triangular, inv

class PreconditionedHessian(LinearOperator):
    def __init__(self, A, L, P):
        super().__init__(shape=L.shape, dtype="double")
        print(L.shape)
        self.L = L
        self.P = P
        self.A = A[:, P]

    def _matvec(self, v):
        print(v.shape)
        v1 = spsolve_triangular(self.L.T, v, lower=False)
        v2 = self.A @ v1
        v3 = self.A.T @ v2
        v4 = spsolve_triangular(self.L, v3, lower=True)
        return v4
