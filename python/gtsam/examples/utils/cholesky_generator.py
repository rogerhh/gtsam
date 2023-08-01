from scikits.sparse.cholmod import cholesky, cholesky_AAt
# from sksparse.cholmod import cholesky, cholesky_AAt
from scipy.sparse import csc_matrix, csr_matrix
from scipy.sparse.linalg import cg, spsolve_triangular, inv
import matplotlib.pyplot as plt
import scipy

class CholeskyGenerator:
    def cholesky(self, Lamb, ridge=0):
        chol_factor = cholesky(Lamb + ridge * scipy.sparse.eye(Lamb.shape[0]))
        L = chol_factor.L()
        P = chol_factor.P()
        return L, P

