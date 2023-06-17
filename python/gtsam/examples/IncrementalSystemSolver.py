import math

import matplotlib.pyplot as plt
import numpy as np

import gtsam
import gtsam.utils.plot as gtsam_plot

from typing import List

from optparse import OptionParser

from scikits.sparse.cholmod import cholesky, cholesky_AAt
# from sksparse.cholmod import cholesky, cholesky_AAt
from scipy.sparse import csc_matrix, csr_matrix
from scipy.sparse.linalg import cg, spsolve_triangular, inv
import matplotlib.pyplot as plt
import scipy

from utils.preconditioner_generator import PreconditionerGenerator
from utils.sparse_linear_system import SparseLinearSystem
from utils.problem_advancer import ProblemAdvancer
from utils.convergence_checker import ConvergenceChecker

from copy import deepcopy

from utils.direct_solver import DirectSolver
from utils.iterative_solver import IterativeSolver
from utils.utils import *

def fillPreconditionerAndPermutation(L, P, A):
    assert(np.allclose(L, np.tril(L)))
    _, old_width = L.shape
    height, width = A.shape
    A = A.toarray()

    P = deepcopy(P)
    P.resize(width)
    P[old_width:width] = np.array(range(old_width, width))

    new_L = np.zeros((width, width))
    new_L[:old_width, :old_width] = L
    for i in range(old_width, width):
        new_L[i, i] = np.linalg.norm(A[i, :])
    
    assert(np.allclose(new_L, np.tril(new_L)))
    return new_L, P

# R should be lower triangular
def applyPreconditioner(A, R):
    A = A.toarray()
    print(R)
    RAT = spsolve_triangular(R, A.T, lower=True)
    RTRAT = spsolve_triangular(R.T, RAT, lower=False)
    ARTR = RTRAT.T
    return ARTR

if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option("-f", "--dataset", dest="dataset",
                      default="", help="Name of dataset")
    parser.add_option("-k", "--K", dest="K",
                      default="1", help="Period of ISAM2 updates")
    parser.add_option("-e", "--epsilon", dest="epsilon",
                      default="0.01", help="Error tolerance")
    parser.add_option("-m", "--max_iter", dest="max_iter",
                      default="10", help="Number of inner loop iterations")
    parser.add_option("-d", "--d_error", dest="d_error",
                      default="0.001", help="If error does not reduce more than d_error, \
                                             we consider it converged")
    parser.add_option("--relinearize_skip", dest="relinearize_skip",
                      default="1", help="Number of steps between relinearization of variable")
    parser.add_option("--print_frequency", dest="print_frequency",
                      default="100", help="Frequency of printing")
    parser.add_option("--num_steps", dest="num_steps",
                      default="100000000", help="Maximum steps")

    parser.add_option("--relin_threshold", dest="relin_threshold",
                      default="0.001", help="Delta norm to relinearize variable")

    (option, args) = parser.parse_args()
    dataset = option.dataset
    K = int(option.K)
    relinearize_skip = int(option.relinearize_skip)
    epsilon = float(option.epsilon)
    d_error = float(option.d_error)
    max_iter = int(option.max_iter)
    print_frequency = int(option.print_frequency)
    num_steps = int(option.num_steps)
    relin_threshold = float(option.relin_threshold)

    dataset_name = gtsam.findExampleDataFile(dataset)
    measurements = gtsam.NonlinearFactorGraph()
    # Priod on the first variable. Add it to factor graph for uniformity
    zero_prior = gtsam.PriorFactorPose2(0, gtsam.Pose2(0, 0, 0), \
                                        gtsam.noiseModel.Unit.Create(3))
    measurements.push_back(zero_prior)
    (dataset_measurements, initial) = gtsam.load2D(dataset_name)
    measurements.push_back(dataset_measurements)

    padv = ProblemAdvancer(measurements)
    spls = SparseLinearSystem(measurements)

    pgen = PreconditionerGenerator(padv, spls, relinearize_skip, relin_threshold)

    problem_type = "iter-refine"
    solver_type = "iterative"

    estimate = gtsam.Values()

    # Only loop through subsampled steps that incremental Cholesky doesn't do well
    for pg_step, end_problem_step in [(800, 804)]:
        # Get optimized problem and solution
        pg_theta, pg_delta, _ = pgen.getThetaAtStep(pg_step, estimate, 1000)
        # Get A and b from optimized solution
        spls.linearizeAll(pg_theta)
        A, b = spls.getSystem()
        Atb = A.T * b
        Lamb = A.T @ A

        chol_factor = cholesky_AAt(A.T)
        L = chol_factor.L()
        P = chol_factor.P()

        theta = pg_theta
        delta = pg_delta
        delta_vec = deepcopy(delta.vector())
        estimate = theta.retract(delta)

        if solver_type == "direct":
            solver = DirectSolver()
        elif solver_type == "iterative":
            solver = IterativeSolver(L, P, tol=1e-10)
        else:
            raise NotImplementedError

        # Get new problem and solution for each new step after the preconditioning step
        for problem_step in range(pg_step + 1, end_problem_step + 1):

            last_measurement_index = spls.measurement_index

            # Get new theta by adding new variables based on current best estimates
            new_nfg, new_theta, new_measurement_index = padv.advanceToStep(problem_step, estimate)
            pgen.update_isam(new_nfg, new_theta)

            theta.insert(new_theta)
            estimate.insert(new_theta)
            spls.addVariables(new_theta)
            spls.addFactors(new_measurement_index)

            convergence_checker = ConvergenceChecker(max_iter=10, atol=0.0017, rtol=1e-4)

            problem_step_iter = -1

            while True:

                problem_step_iter += 1

                # Get new A and b from theta
                spls.linearizeNew(last_measurement_index, new_measurement_index, theta)
                new_A, new_b = spls.getSystem()

                A.resize(new_A.shape)
                A_diff = new_A - A

                new_Lamb = new_A.T @ new_A
                new_Atb = new_A.T @ new_b

                preupdate_error = chi2_red(padv.graph, estimate, padv.factor_dim)
                print("pre update chi2_graph = ", preupdate_error)

                # Solve new system
                if problem_type == "one-step":

                    # Solve AtA delta = Atb 
                    if solver_type == "direct":
                        delta_vec = solver.solve(new_Lamb, new_Atb)

                    elif solver_type == "iterative":
                        new_Atb = new_Atb.toarray()
                        x0 = deepcopy(delta.vector())
                        _, new_A_width = new_A.shape
                        x0.resize((new_A_width, ))
                        delta_vec = solver.solve(new_Lamb, new_Atb, x0=x0)
                    
                elif problem_type == "iter-refine":
                    # Solve A'tA' d = \nu - \Omega x
                    Lamb.resize(new_Lamb.shape)
                    Lamb_diff = new_Lamb - Lamb

                    # DEBUG Start
                    permuted_old_Lamb = L @ L.T
                    permuted_new_Lamb = new_Lamb.toarray()[P[:, np.newaxis], P[np.newaxis, :]]
                    permuted_old_Lamb.resize(permuted_new_Lamb.shape)
                    old_lamb_diff = permuted_new_Lamb - permuted_old_Lamb

                    # print("Rank of new Lambda: ", np.linalg.matrix_rank(new_Lamb.toarray()))
                    # print("Rank of Omega: ", np.linalg.matrix_rank(Lamb_diff.toarray()))
                    # print("Rank of Omega compared to old Lambda: ", np.linalg.matrix_rank(old_lamb_diff))
                    # DEBUG End

                    Atb.resize(new_Atb.shape)
                    Atb_diff = new_Atb - Atb

                    delta_vec = deepcopy(delta_vec)
                    delta_vec.resize(Atb.shape)

                    # need to transform rhs to csc matrix
                    rhs = np.squeeze(Atb_diff.toarray() - (Lamb_diff @ delta_vec))
                    row_indices = range(len(rhs))
                    col_indices = [0 for _ in row_indices]
                    rhs_csc = csc_matrix((rhs, (row_indices, col_indices)))

                    if solver_type == "direct":
                        d = solver.solve(new_Lamb, rhs_csc)
                        delta_vec += d

                        Lamb = new_Lamb
                        Atb = new_Atb

                    elif solver_type == "iterative":

                        # DEBUG Start
                        new_L, new_P = fillPreconditionerAndPermutation(L.toarray(), P, new_A)
                        permuted_new_Lamb = new_Lamb.toarray()[new_P[:, np.newaxis], new_P[np.newaxis, :]]
                        permuted_old_Lamb.resize(permuted_new_Lamb.shape)
                        old_lamb_diff = permuted_new_Lamb - permuted_old_Lamb
                        # s, v, d = np.linalg.svd(np.linalg.inv(new_L @ new_L.T) @ old_lamb_diff)
                        print("Computing svd")
                        s, v0, d = np.linalg.svd(permuted_new_Lamb)
                        s, v1, d = np.linalg.svd(np.linalg.inv(new_L @ new_L.T) @ permuted_new_Lamb)
                        print("Done computing svd")

                        s1 = deepcopy(s[-1])
                        s1[new_P] = s[-1]
                        print(s1)

                        # Get number of distinct singular values
                        v_set = set()
                        for v in v1:
                            v_set.add(round(v, 4))
                        print(f"Number of distinct singular values: {len(v_set)}")

                        fig = plt.figure()
                        plt.plot(range(len(v0)), v0, 'b')
                        plt.plot(range(len(v1)), v1, 'r')
                        plt.yscale("log")
                        plt.title(f"Singular Values Distribution w/ matrix_size={permuted_new_Lamb.shape[0]}")
                        plt.legend(["Unpreconditioned", "Preconditioned"])
                        plt.savefig(f"svd_{problem_step}_{problem_step_iter}.png")
                        plt.show()
                        # AP = new_A.toarray()[:, new_P]
                        # print(new_L.shape, AP.shape)
                        # # print(np.linalg.cond((inv(L @ L.T) @ AP.T @ AP)))
                        # print(np.linalg.cond((new_A.T @ new_A).toarray()))
                        # print(np.linalg.cond(new_A.toarray()))
                        # APL = spsolve_triangular(new_L, AP.T).T
                        # print(np.linalg.cond((APL.T @ APL)))
                        # print(np.linalg.cond((APL)))

                        # DEBUG End

                        # d = solver.solve(new_Lamb, rhs, x0=delta_vec, old_M=Lamb)
                        d = solver.solve(new_Lamb, rhs, x0=None, old_M=Lamb)
                        delta_vec += d

                        Lamb = new_Lamb
                        Atb = new_Atb

                else:
                    raise NotImplementedError

                assert(np.allclose((new_A.T @ new_A @ delta_vec), (new_A.T @ new_b).toarray(), atol=1e-4))
                # Check linearization and update vector values
                relin_factors = set()
                relin_keys = set()
                delta = gtsam.VectorValues()
                theta_update = gtsam.VectorValues()
                var_size = 3
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

                print("relin keys = ", relin_keys)
                print("relin factors = ", relin_factors)
                spls.linearizeSet(relin_factors, theta)

                estimate = theta.retract(delta)
                converged, chi2_graph = convergence_checker.check_convergence(padv.graph, estimate, padv.factor_dim)

                print("post update chi2_graph = ", chi2_graph)

                # input()

                if converged:
                    break
                


