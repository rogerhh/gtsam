"""
problem_advancer.py: Given an initial_estimate, advance the problem until some given step, adding all the measurements to a nonlinear factor graph and all new variables to a vector of values. Does not support retracting measurements. This is because we want to simulate an optimized solution at the preconditioner step
"""

import gtsam
import re

class ProblemAdvancer:
    def __init__(self, measurements):
        self.measurements = measurements
        self.measurement_index = 0
        self.cur_step = 0
        self.factor_dim = 0
        self.graph = gtsam.NonlinearFactorGraph()
        self.key_to_factor_indices = {}

    def add_factor(self, new_nfg, factor):
        factor_index = self.graph.size()
        self.factor_dim += factor.dim()
        new_nfg.push_back(factor)
        self.graph.push_back(factor)

        for key in factor.keys():
            if key not in self.key_to_factor_indices.keys():
                self.key_to_factor_indices[key] = set()

            self.key_to_factor_indices[key].add(factor_index)

    def get_pose2_at(self, estimate, new_theta, step):
        if new_theta.exists(step):
            return new_theta.atPose2(step)
        else:
            assert(estimate.exists(step), f"key = {step}, estimate size = {estimate.size()}")
            return estimate.atPose2(step)

    def advanceToStep(self, end_step, estimate: gtsam.Values):
        new_nfg = gtsam.NonlinearFactorGraph()
        new_theta = gtsam.Values()

        for step in range(self.cur_step , end_step + 1):

            if step == 0:
                new_theta.insert(0, gtsam.Pose2())

            while self.measurement_index < self.measurements.size():
                measurementf = self.measurements.at(self.measurement_index)

                if isinstance(measurementf, gtsam.BetweenFactorPose2):

                    measurement = measurementf
                    key1, key2 = tuple(measurement.keys())

                    if key1 > step or key2 > step:
                        break

                    if key1 != step and key2 != step:
                        measurement.print()
                        raise "Problem in data file, out-of-sequence measurements"

                    # Add a new factor
                    # new_nfg.push_back(factor)
                    self.add_factor(new_nfg, measurement)

                    (cur_key, lower_key) = (key1, key2) if key1 > key2 else (key2, key1)

                    assert(cur_key == step)

                    # Initialize the new variables

                    if not new_theta.exists(step):
                        prev_pose = self.get_pose2_at(estimate, new_theta, lower_key)

                        inverted = (key1 > key2)
                        pose_diff = measurement.measured()
                        if inverted:
                            pose_diff = pose_diff.inverse()

                        new_pose = prev_pose * pose_diff
                        new_theta.insert(step, new_pose)

                        prev_pose = new_pose

                elif isinstance(measurementf, gtsam.PriorFactorPose2):
                    measurement = measurementf
                    self.add_factor(new_nfg, measurement)

                elif isinstance(measurementf, gtsam.BearingRangeFactor2D):
                    measurement = measurementf
                    key1, key2 = tuple(measurement.keys())

                    if key1 > step or key2 > step:
                        break

                    if key1 != step and key2 != step:
                        measurement.print()
                        raise "Problem in data file, out-of-sequence measurements"

                    # Add a new factor
                    # new_nfg.push_back(factor)
                    self.add_factor(new_nfg, measurement)

                    if not new_theta.exists(key2) and key2 == step:
                        pose = self.get_pose2_at(estimate, new_theta, key1)
                        # pose = new_theta.atPose2(key1)
                        measured_bearing = measurement.measured().bearing()
                        measured_range = measurement.measured().range()
                        lm_pose = pose.transformFrom(measured_bearing.rotate(gtsam.Point2(measured_range, 0.0)))
                        new_theta.insert(key2, lm_pose)

                else:
                    print("Unknown factor type: ", type(measurementf))
                    raise NotImplementedError

                self.measurement_index += 1

        self.cur_step = end_step + 1

        return new_nfg, new_theta, self.measurement_index

