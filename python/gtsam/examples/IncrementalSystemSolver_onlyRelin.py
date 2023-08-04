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
from utils.cholesky_generator import CholeskyGenerator

from copy import deepcopy

from utils.preconditioner_update import ExtendedDiagonalPreconditionerUpdater
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
    parser.add_option("--solver", dest="solver_type",
                      default="direct", help="Delta norm to relinearize variable")
    parser.add_option("--ridge", dest="ridge",
                      default="0", help="Delta norm to relinearize variable")
    parser.add_option("--lc_steps_file", dest="lc_steps_file",
                      default="", help="File listing all the loop closure steps.")
    parser.add_option("--lc_lookahead", dest="lc_lookahead",
                      default="1", help="How many steps before the LC step to run factorization.")

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
    ridge = float(option.ridge)
    lc_steps_file = option.lc_steps_file
    lc_steps = readLCSteps(lc_steps_file)
    lc_lookahead = int(option.lc_lookahead)

    dataset_name = gtsam.findExampleDataFile(dataset)
    measurements = gtsam.NonlinearFactorGraph()
    # Prior on the first variable. Add it to factor graph for uniformity
    zero_prior = gtsam.PriorFactorPose2(0, gtsam.Pose2(0, 0, 0), \
                                        gtsam.noiseModel.Unit.Create(3))
    measurements.push_back(zero_prior)
    (dataset_measurements, initial) = gtsam.load2D(dataset_name)
    measurements.push_back(dataset_measurements)

    padv = ProblemAdvancer(measurements)
    spls = SparseLinearSystem(measurements)

    pgen = PreconditionerGenerator(padv, spls, relinearize_skip, relin_threshold)

    cholgen = CholeskyGenerator()

    solver_type = option.solver_type

    estimate = gtsam.Values()

    print(f"dataset: {dataset}\nlc_lookahead: {lc_lookahead}\n")

    # Only loop through subsampled steps that incremental Cholesky doesn't do well
    for lc_step in lc_steps:

        if lc_step > 1000:
            break

        lookahead_step = lc_step - lc_lookahead
        # for pg_step, end_problem_step in [(start_step, start_step + 4)]:

        # Get optimized problem and solution
        lookahead_theta, lookahead_delta, lookahead_estimate = pgen.getThetaAtStep(lookahead_step, estimate, 1000)
        # Get A and b from optimized solution
        spls.linearizeAll(lookahead_theta)
        A_old, b_old = spls.getSystem()
        Atb_old = A_old.T * b_old
        Lamb_old = A_old.T @ A_old

        L_old, P_old = cholgen.cholesky(Lamb_old, 0)

        lc_theta, lc_delta, _ = pgen.getThetaAtStep(lc_step, lookahead_estimate, 1)

        theta = lc_theta
        delta = lc_delta
        delta_vec = deepcopy(delta.vector())
        estimate = theta.retract(delta)
        spls.linearizeAll(lc_theta)
        A_new, b_new = spls.getSystem()
        # Profiling only relin updates
        A_new = A_new[:A_old.shape[0],:A_old.shape[1]]
        b_new = b_new[:b_old.shape[0]]
        # End Profiling only relin updates
        Atb_new = A_new.T * b_new
        Lamb_new = A_new.T * A_new

        print("Matrix size: ", Lamb_new.shape)
        print("Update rank: ", np.linalg.matrix_rank((Lamb_new - Lamb_old).A))

        P_new = augmentPermutation(A_new, P_old)

        pu = ExtendedDiagonalPreconditionerUpdater()
        L_new = pu.updatePreconditioner(L_old, Lamb_old, Lamb_new)

        if solver_type == "direct":
            assert(0)
            solver = DirectSolver()
            new_Lamb = A_new.T @ A_new
            new_Atb = A_new.T @ b_new
            delta_vec = solver.solve(new_Lamb, new_Atb)

        elif solver_type == "iterative":
            A_cond, w_cond = applyPreconditioner(A_new, b_new, L_new, P_new) 
            print("Condition num before conditioning: ", np.linalg.cond(A_new.A))
            print("Condition num after conditioning: ", np.linalg.cond(A_cond))
            # u, d, vt = np.linalg.svd(A_cond)
            # fig = plt.figure()
            # ax = fig.add_subplot(1, 1, 1)
            # ax.set_yscale('log')
            # ax.plot(range(len(d)), d)
            # plt.show()
            Lamb_cond = A_cond.T @ A_cond
            reset_cg_count()
            LP_delta_vec, info = cg(Lamb_cond, w_cond, callback=cg_increment, tol=1e-10)
            print_cg_count()

        else:
            raise NotImplementedError

        continue


        last_measurement_index = spls.measurement_index

        # Get new theta by adding new variables based on current best estimates
        new_nfg, new_theta, new_measurement_index = padv.advanceToStep(lc_step, estimate)
        pgen.update_isam(new_nfg, new_theta)

        theta.insert(new_theta)
        estimate.insert(new_theta)
        spls.addVariables(new_theta)
        spls.addFactors(new_measurement_index)
        spls.linearizeNew(last_measurement_index, new_measurement_index, theta)

        delta_vec.resize(theta.dim())

        convergence_checker = ConvergenceChecker(max_iter=10, atol=0.0017, rtol=1e-4)

        problem_step_iter = -1

        while True:

            problem_step_iter += 1

            # theta, delta should be updated to account for relin
            theta, delta, relin_keys, relin_factors = updateThetaRelin(theta, delta_vec, spls, padv, relin_threshold)
            print("Relin keys: ", relin_keys)
            # Get new A and b from theta
            spls.linearizeSet(relin_factors, theta)
            A_new, b_new = spls.getSystem()

            A_resize, A_tilde, A_prime, A_prime_raw, b_resize, b_tilde, b_prime, b_prime_raw = splitUpdate(A, A_new, b, b_new)
            
            preupdate_error = chi2_red(padv.graph, estimate, padv.factor_dim)
            print("pre update chi2_graph = ", preupdate_error)

            # sqrt of ridge because we're applying to 
            A_new, b_new = applyRidge(A_new, b_new, math.sqrt(ridge))

            # Solve AtA delta = Atb 
            if solver_type == "direct":
                print("condition number of A orig = ", np.linalg.cond(A_new.A))
                new_Lamb = A_new.T @ A_new
                new_Atb = A_new.T @ b_new

                delta_vec = solver.solve(new_Lamb, new_Atb)

            elif solver_type == "iterative":
                P_new = augmentPermutation(A_new, P)
                L_new = generatePreconditioner(A_prime, L, P_new, ridge)
                Lamb_conditioned, w_conditioned = applyPreconditioner(A_new, b_new, L_new, P_new)
                
                LP_delta_vec, info = cg(Lamb_conditioned, w_conditioned, callback=cg_increment, tol=1e-10)

                # global cg_count
                # cg_count = 0
                P_delta_vec = scipy.sparse.linalg.spsolve_triangular(L_new.T, LP_delta_vec, lower=False)
                delta_vec = np.zeros_like(P_delta_vec)
                delta_vec[P_new] = P_delta_vec
            
            # assert(np.allclose((A_new.T @ A_new @ delta_vec), (A_new.T @ b_new).toarray(), atol=1e-3))

            delta = getDelta(delta_vec, spls)
            estimate = theta.retract(delta)

            converged, chi2_graph = convergence_checker.check_convergence(padv.graph, estimate, padv.factor_dim)

            print("post update chi2_graph = ", chi2_graph)
            # plot_poses2d(estimate, plot=True, params={"step": lc_step, "ridge_constant": ridge})

            input()

            if converged:
                break
            


