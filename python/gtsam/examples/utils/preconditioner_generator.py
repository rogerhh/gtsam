"""
Author: rogerhh
preconditioner_generator.py: 
    Generate theta for a solved optimization problem.
    Theta is the linearization point that generates the matrix A and the rhs b
    Delta is the prior solution that solves A^TA delta = A^T b
                             
"""

import gtsam
from utils.utils import *
from utils.problem_advancer import ProblemAdvancer
from utils.sparse_linear_system import SparseLinearSystem

class PreconditionerGenerator:

    def __init__(self, padv: ProblemAdvancer, spls: SparseLinearSystem, 
                       relinearize_skip=1, relin_threshold=0.001):
        self.parameters = gtsam.ISAM2Params()
        self.parameters.setRelinearizeThreshold(relin_threshold)
        self.parameters.setRelinearizeSkip(relinearize_skip)
        self.isam2 = gtsam.ISAM2(self.parameters)

        self.padv = padv
        self.spls = spls
        self.start_step = 0
        self.factor_dim = 0

        self.theta = gtsam.Values()
        self.delta = gtsam.VectorValues()
        self.estimate = gtsam.Values()

        self.nfg = gtsam.NonlinearFactorGraph()

        # self.measurements = measurements
        # self.measurement_index = 0
        # self.theta = gtsam.Values()

    # Add variables until max_step, then optimize until convergence
    def getThetaAtStep(self, max_step, estimate: gtsam.Values, update_period=1000):

        new_nfg = gtsam.NonlinearFactorGraph()
        new_theta = gtsam.Values()

        update_count = 0

        step = self.start_step
        while step < max_step:
            step += update_period
            if step >= max_step:
                step = max_step

            new_nfg, new_theta, max_measurement_index = self.padv.advanceToStep(step, estimate)
            self.factor_dim = self.padv.factor_dim

            self.nfg.push_back(new_nfg)

            self.spls.addVariables(new_theta)
            self.spls.addFactors(max_measurement_index)

            print(f"preconditioner_generator, step = {step}")
            self.isam2.update(new_nfg, new_theta)
            new_nfg = gtsam.NonlinearFactorGraph()
            new_theta = gtsam.Values()

            for _ in range(5):
                estimate = self.isam2.calculateBestEstimate()
                graph = self.isam2.getFactorsUnsafe()
                chi_graph = chi2_red(graph, estimate, self.factor_dim)
                self.isam2.update()

            estimate = self.isam2.calculateBestEstimate()
            delta = self.isam2.getDelta()

            graph = self.isam2.getFactorsUnsafe()
            chi_graph = chi2_red(graph, estimate, self.factor_dim)
            print(chi_graph)


        self.start_step = max_step

        self.theta = self.isam2.getLinearizationPoint()
        self.delta = self.isam2.getDelta()
        self.estimate = estimate

        return self.theta, self.delta, self.estimate

    # Need to be able to update isam even if we use another method to solve
    def update_isam(self, new_nfg, new_theta):
        self.isam2.update(new_nfg, new_theta)

    def get_preconditioner(self):
        pass



