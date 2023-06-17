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

from utils.factor import Factor
from utils.variable import Variable
from copy import deepcopy

# A factor should map to a set of row/column indices, as repeated entries in the csr_matrix constructor gets summed into one entry

def add_factor(factor: gtsam.NonlinearFactor):
    Factor.add_factor(factor)
    # factors.append(factor)
    # nfg.push_back(factor)

def add_variable(key: int, width: int, initial):
    Variable.add_variable(key, width, initial)
    # variables.insert(key, initial)
    # theta.insert(key, initial)

def chi2_red(config: gtsam.Values):
    error = Factor.nfg.error(config)
    graph_dim = len(Factor.b_data)
    dof = graph_dim - config.dim()
    if dof != 0:
        return 2 * error / dof
    else:
        return error

cg_count = 0
def cg_count_increment(xk):
    global cg_count
    cg_count += 1

def conjugate_gradient(A, b, x_0, preconditioner, permutation):
    print(b.shape)
    global cg_count
    factor = cholesky(A)
    R = preconditioner
    P = permutation
    # R = factor.L()
    # P = factor.P()
    PA = A.toarray()[P[:, np.newaxis], P[np.newaxis, :]]
    # print(P)
    assert(np.all(np.linalg.eigvals(PA) > 0))
    # assert(np.allclose((R @ R.T).toarray(), PA))
    Pb = b[P]
    # RA = spsolve_triangular(R, PA)
    # assert(np.allclose(R @ RA, PA))
    # RRA = spsolve_triangular(R.T, RA, lower=False)
    # assert(np.allclose(R.T @ RRA, RA))
    # Rb = spsolve_triangular(R, Pb)
    # assert(np.allclose(R @ Rb, Pb))
    # RRb = spsolve_triangular(R.T, Rb, lower=False)
    # assert(np.allclose(R.T @ RRb, Rb))
    
    # assert(np.allclose(RRA.T, RRA, rtol=1e-20))
    # assert(np.all(np.linalg.eigvals(RRA) > 0))

    cg_count = 0
    maxiter = 10 if len(P) < 50 else len(P) // 10
    x, info = cg(PA, Pb, callback=cg_count_increment, M=inv(R * R.T), maxiter=maxiter, tol=1e-4)
    # x, info = cg(PA, Pb, callback=cg_count_increment, M=inv(R * R.T), tol=1e-4)
    print("cg count = ", cg_count)
    # assert(info == 0)



    # x2 = np.squeeze(np.linalg.inv(RRA) @ RRb)
    # print(x.shape)
    # print(x2.shape)
    # assert(np.allclose(x, x2))
    # RRb = np.squeeze(RRb)
    # x3 = RRA @ x

    # x4 = np.squeeze(PA @ x)
    # print(x4.shape)
    # print(Pb.shape)
    # Pb = np.squeeze(Pb)
    # assert(np.allclose(x4, Pb))

    x5 = deepcopy(x)
    x5[P] = x
    x = x5
    x6 = np.squeeze(A @ x)
    b = np.squeeze(b)

    # for i in range(len(x6)):
    #     print(x6[i], b[i], x6[i] - b[i])

    # assert(np.allclose(A @ x, b, atol=1e-6))
    # L = cholesky(A)
    # x = L.solve_A(b)
    print("cg increment = ", cg_count, " info = ", info)
    # print(x)
    return x

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
    (measurements, initial) = gtsam.load2D(dataset_name)

    # Create iSAM2 parameters which can adjust the threshold necessary to force relinearization and how many
    # update calls are required to perform the relinearization.
    parameters = gtsam.ISAM2Params()
    parameters.setRelinearizeThreshold(0.001)
    parameters.relinearizeSkip = 1
    isam2 = gtsam.ISAM2(parameters)

    measurement_index = 0
    K_count = 0
    prev_pose = gtsam.Pose2()
    # delta = gtsam.VectorValues()    
    # theta = gtsam.Values()                  # theta is the linearization point
    # variables = gtsam.Values()              # variables is theta \oplus delta
    #                                         # and is the current best estimate
    # nfg = gtsam.NonlinearFactorGraph()
    # factors = []                            # Keep a list of factors because it seems like
    #                                         # the python API does not expose random access
    #                                         # to factors
    # key_factor_indices = [[]]

    orderingToKey = []
    keyToOrdering = []

    last_lambda = csc_matrix((0, 0))
    last_Atb = csc_matrix((0, 0))
    last_delta = np.ndarray(shape=(0, 0))
    preconditioner = None
    permutation = None

    for step in range(0, measurements.size()):
        add_variable(step, 3, gtsam.Pose2())

        if step == 0:
            factor = gtsam.PriorFactorPose2(0, gtsam.Pose2(0, 0, 0), \
                                            gtsam.noiseModel.Unit.Create(3))
            add_factor(factor)

        while measurement_index < measurements.size():
            measurementf = measurements.at(measurement_index)

            if isinstance(measurementf, gtsam.BetweenFactorPose2):
                measurement = measurementf
                key1, key2 = tuple(measurement.keys())

                if key1 > step or key2 > step:
                    break

                if key1 != step and key2 != step:
                    measurement.print()
                    raise "Problem in data file, out-of-sequence measurements"

                # Add a new factor
                add_factor(measurement)

                (cur_key, lower_key) = (key1, key2) if key1 > key2 else (key2, key1)

                assert(cur_key == step)

                # Initialize the new variables

                if lower_key == step - 1:
                    inverted = (key1 > key2)
                    pose_diff = measurement.measured()
                    if inverted:
                        pose_diff = pose_diff.inverse()

                    prev_pose = Variable.theta.atPose2(step - 1)
                    new_pose = prev_pose * pose_diff
                    Variable.variables[cur_key].set_initial(new_pose)

            else:
                raise NotImplementedError

            measurement_index += 1

        # Update iSAM2
        if False:
        # if K_count == K or measurement_index == measurements.size():
            K_count = 0;

            # Outerloop nonlinear optimization problem
            nonlin_opt_count = 0
            while True: 
                delta_copy = gtsam.VectorValues()
                for variable in Variable.variables:
                    delta_norm = variable.delta_norm()

                    if delta_norm >= relin_threshold:
                        delta_copy.insert(variable.key, Variable.delta.at(variable.key))
                        for factor in variable.factors:
                            Factor.relin_factors.add(factor)

                Variable.theta = Variable.theta.retract(delta_copy)

                if len(Factor.relin_factors) == 0:
                    break
                else: 
                    nonlin_opt_count += 1

                for factor in Factor.relin_factors:
                    # print("Relinearize factor ", factor.factor_index)
                    factor.linearize(Variable.theta)
                Factor.relin_factors.clear()

                A_matrix = csr_matrix((Factor.A_data, (Factor.A_rows, Factor.A_cols)))
                b_matrix = csr_matrix((Factor.b_data, (Factor.b_rows, Factor.b_cols)))

                lambda_m = A_matrix.T * A_matrix
                Atb = A_matrix.T * b_matrix

                last_lambda.resize(lambda_m.shape)
                lambda_diff = lambda_m - last_lambda
                last_lambda = lambda_m
                last_Atb.resize(Atb.shape)
                Atb_diff = Atb - last_Atb
                last_Atb = Atb

                # print("lambda_diff = ", lambda_diff)
                # print("Atb_diff = ", Atb_diff)

                if step % 10 == 0 or preconditioner is None:
                    L = cholesky(lambda_m)
                    delta = L.solve_A(Atb.toarray())
                    # L = cholesky_AAt(A_matrix.T)
                    # delta = L.solve_A(Atb)
                    last_delta.resize(delta.shape)
                    delta_diff = last_delta - delta
                    last_delta = delta
                    # print("delta_diff = ", delta_diff)
                    Variable.update_delta(delta)
                    preconditioner = L.L()
                    permutation = deepcopy(L.P())

                else:
                    old_len_rows = preconditioner.shape[0]
                    new_len_rows = lambda_m.shape[0]
                    diag = lambda_m.diagonal()
                    diag_indices = range(old_len_rows, new_len_rows)
                    diag_root = diag[old_len_rows:new_len_rows] ** 0.5
                    diag_root = [(x if x >= 1 else 1) for x in diag_root]
                    diag_precond = csc_matrix((diag_root, (diag_indices, diag_indices)), \
                                               shape=lambda_m.shape)
                    preconditioner.resize(lambda_m.shape)
                    preconditioner = preconditioner + diag_precond
                    permutation.resize(new_len_rows)
                    permutation[old_len_rows:new_len_rows] = np.array(range(old_len_rows, new_len_rows))
                    last_delta = np.resize(last_delta, Atb.shape)
                    d = conjugate_gradient(lambda_m, Atb_diff.toarray() 
                                                     - lambda_diff.toarray() @ last_delta, 
                                           csc_matrix(Atb.shape),
                                           preconditioner,
                                           permutation)
                    print(last_delta.shape)
                    print(d.shape)
                    delta = np.squeeze(last_delta) + d
                    last_delta = delta
                    print(type(delta))
                    print(delta.shape)
                    Variable.update_delta(np.squeeze(delta))

                error = chi2_red(Variable.get_best_estimate())
                print(f"step = {step}, error = {error}")

                if error < 0.0018:
                    break

            print("Number of nonlinear optimization steps = ", nonlin_opt_count)

        K_count += 1

        if step % print_frequency == 0:
            # error = nfg.error(variables)
            error = Factor.error(Variable.theta)
            print(f"step = {step}, error = {error}")
            Variable.print_variables()

        if step >= num_steps:
            break


    error = chi2_red(Variable.get_best_estimate())
    print(f"Final error = {error}")
    Variable.print_variables()

