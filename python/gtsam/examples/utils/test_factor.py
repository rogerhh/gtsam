from factor import Factor
from variable import Variable
import gtsam
import numpy as np

if __name__ == "__main__":
    p1 = gtsam.Pose2(0, 0, 0)
    p2 = gtsam.Pose2(1, 1, np.pi)

    pd = p1.between(p2)
    
    Variable.add_variable(0, 3, p1)
    Variable.add_variable(1, 3, p2)

    UNIT_NOISE = gtsam.noiseModel.Unit.Create(3)

    gtsam_factor = gtsam.BetweenFactorPose2(0, 1, pd, UNIT_NOISE)

    Factor.add_factor(gtsam_factor)

    print(Factor.A_max_len)
    print(Factor.b_max_len)
    factor = Factor(gtsam_factor)

    print(Factor.A_max_len)
    print(Factor.b_max_len)
