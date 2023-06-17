"""
utils.py: Utility functions that may be used by all files
Author: rogerhh
Date: 2023-05-30
"""

import typing

import gtsam

def chi2_red(nfg: gtsam.NonlinearFactorGraph, config: gtsam.Values, factor_dim) -> float:
    error = nfg.error(config)
    graph_dim = factor_dim
    dof = graph_dim - config.dim()
    assert(dof >= 0)
    if dof != 0:
        return 2 * error / dof
    else:
        print("dof == 0!")
        raise RuntimeError


cg_count = 0
def cg_increment(x):
    global cg_count
    cg_count += 1
    print(cg_count)
