"""
utils.py: Utility functions that may be used by all files
Author: rogerhh
Date: 2023-05-30
"""

import typing

import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
# from scikits.sparse.cholmod import cholesky, cholesky_AAt
from sksparse.cholmod import cholesky, cholesky_AAt
from scipy.sparse import csc_matrix, csr_matrix
from scipy.sparse.linalg import cg, spsolve_triangular, inv
import scipy

import gtsam

def chi2_red(nfg: gtsam.NonlinearFactorGraph, config: gtsam.Values, factor_dim) -> float:
    error = nfg.error(config)
    graph_dim = factor_dim
    dof = graph_dim - config.dim()
    assert(dof >= 0)
    if dof != 0:
        return 2 * error / dof
    else:
        print("dof == 0!")
        raise RuntimeError


cg_count = 0
def cg_increment(x):
    global cg_count
    cg_count += 1
    print(cg_count)

def plot_poses2d(poses, filename=None, params={}, plot=True, save=False):
    # plt.figure()
    x = []
    y = []
    for i in range(poses.size()):
        pose = poses.atPose2(i)
        x.append(pose.x())
        y.append(pose.y())
    title = ""
    for key, item in params.items():
        title += f"{key}={item} "
    plt.title(title)
    plt.plot(x, y, '-', alpha=0.5, color="green")
    if plot:
        plt.show()
    if save:
        asserT(filename is not None)
        plt.savefig(filename)


# Update delta and theta based on relin threshold
def updateThetaRelin(theta, delta_vec, spls, padv, relin_threshold=0.001):

    # Check linearization and update vector values
    relin_factors = set()
    relin_keys = set()
    delta = gtsam.VectorValues()
    theta_update = gtsam.VectorValues()
    for key, (col, width) in enumerate(spls.key_to_col):
        delta_key = delta_vec[col:col+width]
        delta_norm = np.linalg.norm(delta_key, ord=np.inf)
        if delta_norm >= relin_threshold:
            # This delta will be retracted from theta
            theta_update.insert(key, delta_key)
            delta.insert(key, np.zeros_like(delta_key))
            relin_factors = relin_factors | padv.key_to_factor_indices[key]
            relin_keys.add(key)

        else:
            # This delta will not be retracted from theta
            delta.insert(key, delta_key)

    theta = theta.retract(theta_update)

    return theta, delta, relin_keys, relin_factors

# Update delta and theta without considering relin
def getDelta(delta_vec, spls):

    # Check linearization and update vector values
    relin_factors = set()
    relin_keys = set()
    delta = gtsam.VectorValues()
    for key, (col, width) in enumerate(spls.key_to_col):
        delta_key = delta_vec[col:col+width]
        delta.insert(key, delta_key)

    return delta

def splitUpdate(oldA, newA, oldb, newb):
    h1, w1 = oldA.shape
    h2, w2 = newA.shape

    A_copy = deepcopy(oldA)
    A_copy.resize(newA.shape)
    A_prime = csr_matrix(newA.shape)
    A_prime[h1:h2] = newA[h1:h2]
    A_prime_raw = newA[h1:h2]
    A_tilde = newA - A_prime - A_copy

    assert(oldb.shape[0] == h1)
    assert(newb.shape[0] == h2)

    b_copy = deepcopy(oldb)
    b_copy.resize(newb.shape)
    b_prime = csr_matrix(newb.shape)
    b_prime[h1:h2] = newb[h1:h2]
    b_prime_raw = newb[h1:h2]
    b_tilde = newb - b_prime - b_copy

    return A_copy, A_tilde, A_prime, A_prime_raw, b_copy, b_tilde, b_prime, b_prime_raw

def applyRidge(A, b, ridge_constant=0):
    assert(scipy.sparse.issparse(A))

    A_copy = deepcopy(A)
    height, width = A_copy.shape
    A_copy.resize(height + width, width)
    A_copy.tocsr()[height:] += ridge_constant * scipy.sparse.eye(width)

    b_copy = deepcopy(b)
    b_copy.resize((height + width, 1))

    return A_copy, b_copy

def augmentPermutation(A, P):
    width = A.shape[1]
    old_width = P.shape[0]
    P_new = deepcopy(P)
    P_new.resize((width, ))
    P_new[old_width:width] = np.arange(old_width, width)
    return P_new

def generatePreconditioner(A_prime, L, P, ridge):
    old_width = L.shape[1]
    new_width = A_prime.shape[1]

    print("P = ", P)
    # First apply ridge constant to L
    Lambda_old = L @ L.T
    Lambda_old += ridge * scipy.sparse.eye(old_width)
    chol_factor = cholesky(Lambda_old, ordering_method="natural")
    L = chol_factor.L()
    L.resize((new_width, new_width))
    print("P1 = ", chol_factor.P())

    # Compute Lambda prime with ridge only to the new variables
    A_prime_permuted = A_prime.tocsc()[:, P]
    Lambda_prime = A_prime_permuted.T @ A_prime_permuted
    Lambda_prime[old_width:new_width, old_width:new_width] += ridge * scipy.sparse.eye(new_width - old_width)
    Lambda_prime_diag = Lambda_prime.diagonal()

    L_diag = L.diagonal()
    new_diag = np.sqrt(L_diag * L_diag + Lambda_prime_diag)

    # L = scipy.sparse.eye(new_width)
    # L[range(old_width,new_width), range(old_width,new_width)] = new_diag[old_width:new_width]
    L[range(0,new_width), range(0,new_width)] = new_diag[0:new_width]

    return L

def applyPreconditioner(A, b, L, P):
    A_permuted = A.tocsc()[:, P].A
    print(A_permuted.shape)
    print(L.shape)
    A_conditioned = scipy.sparse.linalg.spsolve_triangular(L, A_permuted.T, lower=True).T
    Lamb_conditioned = A_conditioned.T @ A_conditioned

    w_conditioned = A_conditioned.T @ b

    print("orig A cond = ", np.linalg.cond(A_permuted))
    print("conditioned A cond = ", np.linalg.cond(A_conditioned))

    return Lamb_conditioned, w_conditioned


