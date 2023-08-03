import math

import matplotlib.pyplot as plt
import numpy as np

import gtsam
import gtsam.utils.plot as gtsam_plot

from typing import List

import math

if __name__ == "__main__":

    p1 = gtsam.Pose2(0, 0, 0)
    p2 = gtsam.Pose2(1, 1, np.pi * 2 + 1)
    pd = p1.between(p2)
    pd2 = p1.inverse() * p2
    print("pd2 = ", pd2)

    pose_diff = gtsam.Pose2(0, 0, 0)
    
    theta = gtsam.Values()
    theta.insert(0, p1)
    theta.insert(1, p2)
    factor = gtsam.BetweenFactorPose2(0, 1, pose_diff, gtsam.noiseModel.Unit.Create(3))

    error = factor.error(theta)
    print(error)

    print(p1.matrix())
    print(p2.matrix())
    print(p2.translation())
    print(p2.rotation())
    print(p2.x(), p2.y(), p2.theta())
    print(pd)
    print(p1 * pd)
    print(gtsam.Pose2.Logmap(pd))
    print((np.linalg.norm(gtsam.Pose2.Logmap(pd)) ** 2) / 2)

