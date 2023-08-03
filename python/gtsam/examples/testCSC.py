from scikits.sparse.cholmod import cholesky
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import scipy

if __name__ == "__main__":
    data = [1, 1, 1, 1, 1, 1, 2, 3]
    rows = [0, 1, 2, 3, 4, 5, 3, 4]
    cols = [0, 1, 2, 3, 4, 5, 3, 4]
    matrix = csr_matrix((data, (rows, cols)), shape=(6, 6))
    print(matrix.toarray())
