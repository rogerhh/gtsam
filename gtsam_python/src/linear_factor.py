"""
Author: Roger Hsiao
linear_factor.py: Factor with a linear constraint on the variables. I.e. a Gaussian distribution

"""
from abc import ABC, abstractmethod
from typing import List, Tuple
import numpy as np
import numpy.typing as npt

from gtsam_python.src.factor import Factor

class LinearFactor(Factor):

    def __init__(self, varIndices: List[int], A: List[npt.NDArray], b: npt.NDArray):
        assert(len(varIndices) == len(A))
        super().__init__(varIndices)
        self.A = A
        self.b = b

    # Given an X, linearize this factor to a A matrix and a b vector
    def linearize(self, theta: List[npt.NDArray]) -> Tuple[List[npt.NDArray], npt.NDArray]:
        return self.A, self.b

    # Given an X, calculate the prediction error of this factor
    def error(self, theta: List[npt.NDArray]) -> float:
        y = np.zeros_like(self.b)
        for i, index in enumerate(self.vars):
            assert(index < len(theta))
            A_i = self.A[i]
            theta_i = theta[i]
            assert(A_i.shape[0] == self.b.shape[0])
            assert(A_i.shape[1] == theta_i.shape[0])
            assert(theta_i.shape[1] == 1)
            y += A_i @ theta_i

        norm = np.linalg.norm(y - self.b)
        norm = norm ** 2

        return norm
        
