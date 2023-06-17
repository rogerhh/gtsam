"""
    Author: rogerhh
    Date: May 22, 2023
    factor.py: 
    - A Factor should own a set of indices in the csr_matrix contructor list and a set of Variable's. Through the starting_col and the width of the Variabe's, the Factor owns an exclusive subset of the csr_matrix construction lists
    - Repeated entries in the csr_matrix constructor gets summed into one entry. 
    - Adding new Factor's appends the (row, column) indices and the entries to CSR matrix
    - Relinearizing a Factor changes the underlying values of the entries, and should provide the difference between the new and old values
    - Reordering the variables changes the (row, column) indices the Factor owns

"""

from typing import List
import gtsam
from utils.variable import Variable

class Factor:
    factors = []
    relin_factors = set()
    nfg = gtsam.NonlinearFactorGraph()

    max_rows = 0

    A_max_len = 0
    A_data = []
    A_rows = []
    A_cols = []

    b_max_len = 0
    b_data = []
    b_rows = []
    b_cols = []

    @staticmethod
    def add_factor(factor: gtsam.NonlinearFactor):
        Factor.factors.append(Factor(factor))
        Factor.relin_factors.add(Factor.factors[-1])
        Factor.nfg.push_back(factor)

    @staticmethod
    def error(theta: gtsam.Values):
        return Factor.nfg.error(theta)

    def __init__(self, nlf: gtsam.NonlinearFactor):
        self.nlf = nlf
        self.factor_index = len(Factor.factors)
        self.dim = self.nlf.dim()
        self.start_row = Factor.max_rows
        Factor.max_rows += self.dim

        self.A_start_index = len(Factor.A_data)
        self.A_width = 0
        for key in self.nlf.keys():
            self.A_width += Variable.variables[key].width
        Factor.A_max_len += self.A_width * self.dim
        
        self.b_start_index = len(Factor.b_data)
        Factor.b_max_len += self.dim

        for key in self.nlf.keys():
            Variable.variables[key].factors.append(self)


        col_indices = []

        for key in self.nlf.keys():
            (col, width) = Variable.get_col_and_width(key)
            for j in range(width):
                c = col + j
                for i in range(self.dim):
                    r = self.start_row + i
                    Factor.A_data.append(0)
                    Factor.A_cols.append(c)
                    Factor.A_rows.append(r)
        
        for i in range(self.dim):
            Factor.b_data.append(0)
            Factor.b_rows.append(self.start_row + i)
            Factor.b_cols.append(0)

    def linearize(self, theta: gtsam.Values):
        jacobian_factor = self.nlf.linearize(theta)
        matrixA = jacobian_factor.getA()
        vectorb = jacobian_factor.getb()

        A_index = self.A_start_index
        col_start = 0

        for key in self.nlf.keys():
            (col, width) = Variable.get_col_and_width(key)

            for j in range(width):
                c = col + j
                for i in range(self.dim):
                    r = self.start_row + i

                    assert(Factor.A_rows[A_index] == r)
                    assert(Factor.A_cols[A_index] == c)

                    Factor.A_data[A_index] = matrixA[i, col_start + j]

                    A_index += 1

            col_start += width

        for i in range(self.dim):
            r = self.start_row + i

            Factor.b_data[r] = vectorb[i]


        
        return None
        

