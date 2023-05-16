from gtsam_python.src.linear_factor import LinearFactor
import numpy as np

if __name__ == "__main__":
    A00 = np.array(np.arange(1, 10)).reshape(3, 3)
    A01 = np.array(np.arange(2, 11)).reshape(3, 3)
    b0 = np.array(np.arange(0, 3)).reshape(3, 1)


    factor1 = LinearFactor(varIndices=[0, 1], A=[A00, A01], b=b0)

    x1 = np.array([3, 3, 3]).reshape(3, 1)
    x2 = np.array([4, 4, 4]).reshape(3, 1)
