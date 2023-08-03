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
from utils.iterative_solver_simple import IterativeSolver
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
        new_L[i, i] = np.linalg.norm(A[:, i])

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

# R should be lower triangular
def applyPreconditionerDense(A, R):
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

    problem_type = "one-step"
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

        ridge = 0
        # chol_factor = cholesky_AAt(A.T)
        # L = chol_factor.L()
        # P = chol_factor.P()
        chol_factor = cholesky(Lamb + ridge * scipy.sparse.eye(Lamb.shape[0]))
        L = chol_factor.L()
        P = chol_factor.P()

        theta = pg_theta
        delta = pg_delta
        delta_vec = deepcopy(delta.vector())
        estimate = theta.retract(delta)

        if solver_type == "direct":
            solver = DirectSolver()
        elif solver_type == "iterative":
            default_L = scipy.sparse.eye(Lamb.shape[0])
            default_P = np.array(np.arange(P.shape[0]))
            solver = IterativeSolver(default_L, default_P, tol=1e-10)
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
                        print("condition number orig = ", np.linalg.cond(new_Lamb.toarray()))
                        print("condition number ridge = ", np.linalg.cond((new_Lamb + ridge * scipy.sparse.eye(new_Lamb.shape[0]).toarray())))
                        # delta_vec = solver.solve(new_Lamb, new_Atb)
                        delta_vec = solver.solve(new_Lamb + ridge * scipy.sparse.eye(new_Lamb.shape[0]), new_Atb)

                    elif solver_type == "iterative":
                        new_L, new_P = fillPreconditionerAndPermutation(L.toarray(), P, new_A)
                        permuted_new_Lamb = new_Lamb.toarray()[new_P[:, np.newaxis], new_P[np.newaxis, :]]
                        conditioned_Lamb = applyPreconditionerDense(permuted_new_Lamb, new_L)

                        # print("condition number orig = ", np.linalg.cond(new_Lamb.toarray()))
                        # print("condition number permuted = ", np.linalg.cond(permuted_new_Lamb))
                        # print("condition number conditioned = ", np.linalg.cond(conditioned_Lamb))

                        permuted_A = new_A.toarray()[:, new_P]
                        conditioned_A = scipy.linalg.solve_triangular(new_L, permuted_A.T, lower=True).T
                        assert(np.allclose(conditioned_A @ new_L.T, permuted_A))
                        A_Lamb = permuted_A.T @ permuted_A
                        cond_A_Lamb = conditioned_A.T @ conditioned_A

                        # fig = plt.figure()
                        # ax = fig.add_subplot(1, 1, 1)
                        # ax.set_yscale('log')
                        s, v1, dt = np.linalg.svd(permuted_A)
                        s, v2, dt = np.linalg.svd(conditioned_A)
                        perturb_vec = (dt.T)[:, -1].reshape((-1, 1))
                        perturb_vec = scipy.linalg.solve_triangular(new_L.T, perturb_vec, lower=False)
                        # print(v1[-1], perturb_vec.shape)
                        print(v1[-1], perturb_vec.shape)

                        # plt.plot(np.arange(v1.shape[0]), v1)
                        # plt.plot(np.arange(v2.shape[0]), v2)
                        # plt.show()
                        # exit()

                        print("condition number orig A = ", np.linalg.cond(permuted_A))
                        print("condition number conditioned A = ", np.linalg.cond(conditioned_A))
                        print("condition number A to Lamb = ", np.linalg.cond(A_Lamb))
                        print("condition number conditioned A to Lamb = ", np.linalg.cond(cond_A_Lamb))
                        new_Atb = new_Atb.toarray()
                        x0 = deepcopy(delta.vector())
                        _, new_A_width = new_A.shape
                        x0.resize((new_A_width, ))
                        Py = new_Atb[new_P]
                        # print("Solving original")
                        # P_delta_vec = solver.solve(A_Lamb, Py, x0=x0)
                        # print("Solving conditioned no initial")
                        # P_delta_vec = solver.solve(cond_A_Lamb, Py)
                        print("Solving conditioned")
                        LPy = scipy.linalg.solve_triangular(new_L, Py, lower=True)
                        # LP_delta_vec = solver.solve(cond_A_Lamb, LPy, x0=x0)
                        LP_delta_vec, info = cg(cond_A_Lamb, LPy, callback=cg_increment)
                        P_delta_vec = scipy.linalg.solve_triangular(new_L.T, LP_delta_vec, lower=False)
                        # print("Solving conditioned 2")
                        # conda_A_Lamb = scipy.linalg.solve_triangular(new_L, A_Lamb, lower=True)
                        # conda_A_Lamb = scipy.linalg.solve_triangular(new_L.T, cond_A_Lamb, lower=False).T
                        # LPy = scipy.linalg.solve_triangular(new_L, Py, lower=True)
                        # LP_delta_vec = solver.solve(cond_A_Lamb, LPy, x0=x0)
                        # P_delta_vec = scipy.linalg.solve_triangular(new_L.T, LP_delta_vec, lower=False)
                        scale = 10
                        P_delta_vec += scale * perturb_vec
                        delta_vec = np.zeros_like(P_delta_vec)
                        delta_vec[new_P] = P_delta_vec
                        
                    
                else:
                    raise NotImplementedError

                # assert(np.allclose((new_A.T @ new_A @ delta_vec), (new_A.T @ new_b).toarray(), atol=1e-3))
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
                # print("relin factors = ", relin_factors)
                spls.linearizeSet(relin_factors, theta)

                estimate = theta.retract(delta)
                converged, chi2_graph = convergence_checker.check_convergence(padv.graph, estimate, padv.factor_dim)

                print("post update chi2_graph = ", chi2_graph)

                # plot_poses2d(estimate, f"output/poses_{pg_step}_{scale}.png", step=pg_step, error=chi2_graph, scale=scale)
                print("Done")

                input()

                if converged:
                    break
                


