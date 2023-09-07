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
from utils.sparse_linear_system3D import SparseLinearSystem
from utils.problem_advancer3D import ProblemAdvancer
from utils.convergence_checker import ConvergenceChecker
from utils.cholesky_generator import CholeskyGenerator

from copy import deepcopy

from utils.preconditioner_update import *
from utils.direct_solver import DirectSolver
from utils.iterative_solver import IterativeSolver
from utils.linear_operator import PreconditionedHessian
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
    parser.add_option("--preconditioner_type", dest="pu_type",
                      default="identity", help="What type of preconditioner updater to use.")

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
    print(dataset, dataset_name)
    assert(dataset_name is not None)
    measurements = gtsam.NonlinearFactorGraph()
    # Prior on the first variable. Add it to factor graph for uniformity
    zero_prior = gtsam.PriorFactorPose3(0, gtsam.Pose3(), \
                                        gtsam.noiseModel.Unit.Create(6))
    measurements.push_back(zero_prior)
    (dataset_measurements, initial) = gtsam.load3D(dataset_name)
    measurements.push_back(dataset_measurements)

    padv = ProblemAdvancer(measurements)
    spls = SparseLinearSystem(measurements)

    pgen = PreconditionerGenerator(padv, spls, relinearize_skip, relin_threshold)

    cholgen = CholeskyGenerator()

    solver_type = option.solver_type

    pu_type = option.pu_type
    if pu_type == "trueidentity":
        pu = TrueIdentityPreconditionerUpdater()
    elif pu_type == "identity":
        pu = IdentityPreconditionerUpdater()
    elif pu_type == "extdiag":
        pu = ExtendedDiagonalPreconditionerUpdater()
    elif pu_type == "incomplcholdiag":
        pu = IncompleteCholeskyOnDiagPreconditionerUpdater()
    elif pu_type == "cholnewblock":
        pu = CholeskyOnNewBlockPreconditionerUpdater()
    elif pu_type == "incompletechol":
        pu = IncompleteCholeskyPreconditionerUpdater()
    elif pu_type == "incompletecholrelin":
        pu = IncompleteCholeskyWithRelinLambdaPreconditionerUpdater()
    elif pu_type == "cholupdate":
        pu = CholeskyUpdatePreconditionerUpdater()
    elif pu_type == "selcholupdate":
        pu = SelectiveCholeskyUpdatePreconditionerUpdater()
    elif pu_type == "selcholupdate2":
        pu = SelectiveCholeskyUpdatePreconditionerUpdater2()
    else:
        raise NotImplementedError


    estimate = gtsam.Values()

    print(f"dataset: {dataset}\npreconditioner type: {pu_type}\nlc_lookahead: {lc_lookahead}\n")

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
        Atb_new = A_new.T * b_new
        Lamb_new = A_new.T * A_new

        height_old, width_old = A_old.shape
        height_new, width_new = A_new.shape
        A_prime = A_new[height_old:height_new,:]
        A_tilde = A_new[:height_old,:width_old] - A_old


        print("Old A shape: ", A_old.shape)
        print("New A shape: ", A_new.shape)
        print("Matrix size: ", Lamb_new.shape)
        print("Relin rank: ", np.linalg.matrix_rank((A_tilde.T @ A_tilde).A))
        print("Observe rank: ", np.linalg.matrix_rank(A_prime.A))

        P_new = augmentPermutation(A_new, P_old)

        L_new = pu.updatePreconditioner(L_old=L_old, P_old=P_old, P_new=P_new, \
                                        Lamb_old=Lamb_old, Lamb_new=Lamb_new, \
                                        A_old=A_old, A_prime=A_prime, A_tilde=A_tilde)

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

            linOps = PreconditionedHessian(A_new, L_new, P_new)

            reset_cg_count()
            LP_delta_vec, info = cg(linOps, w_cond, callback=cg_increment, tol=1e-10, maxiter=height_new)
            print_cg_count()

            input()

        else:
            raise NotImplementedError

        continue

