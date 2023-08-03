"""
    Author: rogerhh
    Date: May 22, 2023
    Variable.py: 
    - A Variable should own a start_col and the width of the variable.
    - A Variable should also own its own linearization point and its current best estimation (This can be skipped since we want to have a unified access to theta and delta)
    - Reordering variables changes the starting col member

"""

from typing import List
import numpy as np
import gtsam

class Variable:
    max_col = 0
    variables = []
    theta = gtsam.Values()
    delta = gtsam.VectorValues()
    orderingToKey = []
    keyToOrdering = []

    @staticmethod
    def add_variable(key: int, width: int, initial):
        print(f"add_variable {key}, {width}")
        Variable.variables.append(Variable(key, Variable.max_col, width))
        Variable.max_col += width

        Variable.theta.insert(key, initial)
        Variable.delta.insert(key, np.zeros((width, )))

        Variable.orderingToKey.append(len(Variable.orderingToKey))
        Variable.keyToOrdering.append(len(Variable.keyToOrdering))

    @staticmethod
    def at(key: int):
        return theta[key]

    @staticmethod
    def print_variables():
        Variable.theta.print()

    def __init__(self, key: int, start_col: int, width: int):
        self.key = key
        self.start_col = start_col
        self.width = width
        self.factors = []

    def set_initial(self, initial):
        Variable.theta.update(self.key, initial)

    @staticmethod
    def get_col_and_width(key: int):
        return Variable.variables[key].start_col, Variable.variables[key].width

    @staticmethod
    def ordering_less(key1: int, key2: int) -> bool:
        return Variable.keyToOrdering(key1) < Variable.keyToOrdering(key2)

    @staticmethod
    def update_delta(delta: np.ndarray):
        # Need to make a new delta then copy it over because 
        # direct update is not exposed
        delta_new = gtsam.VectorValues()

        start_row = 0
        for key, order in enumerate(Variable.keyToOrdering):
            var = Variable.variables[key]
            delta_new.insert(key, delta[start_row:start_row + var.width])

            start_row += var.width

        Variable.delta.update(delta_new)

    def update_theta(self):
        print(Variable.theta.exists(self.key))
        print(Variable.theta.atVector(self.key))
        exit()
        new_theta = Variable.theta.at(self.key).retract(Variable.delta.at(self.key))
        Variable.theta.update(self.key, new_theta)
        print(new_theta)
        exit()

    @staticmethod
    def get_best_estimate():
        return Variable.theta.retract(Variable.delta)

    def delta_norm(self):
        return np.linalg.norm(Variable.delta.at(self.key), ord=np.inf)
