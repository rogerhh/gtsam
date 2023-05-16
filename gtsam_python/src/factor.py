"""
Author: Roger Hsiao
factor.py: Defines an abstract base class Factor that denotes a relationship between one or more variables

"""

from abc import ABC, abstractmethod
import numpy as np
import numpy.typing as npt
import typing

class Factor(ABC):

    def __init__(self, varIndices: typing.List[int]):
        self.vars = varIndices

    # Given an X, linearize this factor to a A matrix and a b vector
    @abstractmethod
    def linearize(self, theta: typing.List[npt.NDArray]):
        ...

    # Given an X, calculate the prediction error of this factor
    @abstractmethod
    def error(self, theta: typing.List[npt.NDArray]):
        ...
