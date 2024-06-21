"""
generate_dataset_vio.py
Generate a full dataset for VIO
This means 

Usage: generate_full_tree.py --infile <tree_structure_file> --seed 0 --outfile <outfile_path>
"""

import numpy as np
import scipy
from scipy.stats import special_ortho_group
import sys
from optparse import OptionParser
import math
from copy import deepcopy
from utils import *
import glob
import os
from timestep import Timestep

def read_until(fin, s):
    while True:
        line = fin.readline()
        if not line:
            print(f"Pattern {s} not found!")
            exit(0)

        if s in line:
            break

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("--indir", dest="indir", 
                      default=None, help="The path to the input directory")
    parser.add_option("--outdir", dest="outdir", 
                      default=None, help="The path to the output directory")
    parser.add_option("--num_threads_file", dest="num_threads_file", 
                      default=None, help="The file path to the threads file that we use")
    parser.add_option("--start_step", dest="start_step", type="int",
                      default=2, help="Starting step, inclusive")
    parser.add_option("--end_step", dest="end_step", type="int",
                      default=30, help="Ending step, inclusive")
    parser.add_option("--period", dest="period", type="int",
                      default=1, help="Step size of incrementing timesteps")
    parser.add_option("--no_values", dest="no_values", 
                      action="store_true", help="If the generated dataset has no matrix values")
    parser.add_option("--vio", dest="vio", action="store_true",
                      default=None, help="Set if generating vio dataset. That means, generate dataset without reordering and discard all factors that loop back more than vio_lag steps")
    parser.add_option("--vio_lag", dest="vio_lag", type="int",
                      default=1e10, help="See vio")
    parser.add_option("--vio_scale", dest="vio_scale", type="float",
                      default=1.0, help="Multiple marginal_factor by vio_scale for numerical stability")
    parser.add_option("--relin_thresh", dest="relin_thresh", type="float",
                      default=0.1, help="Relinearize threshold")
    parser.add_option("--seed", dest="seed", type="int",
                      default=None, help="random seed for numpy sampling")
    (options, args) = parser.parse_args()

    indir = options.indir
    outdir = options.outdir
    num_threads_file = options.num_threads_file
    start_step = options.start_step
    end_step = options.end_step
    period = options.period
    no_values = options.no_values
    vio = options.vio
    vio_lag = options.vio_lag
    relin_thresh = options.relin_thresh

    Timestep.vio_scale = options.vio_scale
    Timestep.no_values = no_values

    if options.seed is not None:
        np.random.seed(options.seed)

    step_num_threads = None
    if num_threads_file is not None:
        step_num_threads = []
        with open(num_threads_file, "r") as fin:
            while True:
                line = fin.readline()
                if not line:
                    break
                step_num_threads.append(int(line))
        Timestep.step_num_threads = step_num_threads

        Timestep.max_num_threads = max(step_num_threads)
        Timestep.min_num_threads = min(step_num_threads)


    files = glob.glob(os.path.join(indir, "step-*.out"))

    timesteps = [None for _ in range(0, end_step+1)]
    
    for step in range(start_step, end_step+1):
        if step % period != 0:
            continue

        infile_pattern = f"step-{step}.out"
        infile = None
        pred_infile_pattern = f"step-{step}-pred_cycles.out"
        pred_infile = None

        for filename in files:
            if infile_pattern in filename:
                infile = filename
                break

        for filename in files:
            if pred_infile_pattern in filename:
                pred_infile = filename
                break

        num_threads = -1
        relin_cost = -1
        if pred_infile is not None:
            with open(pred_infile, "r") as fin:
                read_until(fin, "num threads")
                num_threads = int(fin.readline())
                read_until(fin, "relin cost")
                relin_cost = int(fin.readline())

        if infile is None:
            print(f"file at step {step} not found!", step)
            assert(0)
    
        with open(infile, "r") as fin:
            timesteps[step] = Timestep(fin, step, num_threads=num_threads, relin_cost=relin_cost)

    # If vio, remove all factors that span over vio_lag steps
    # Also remove all factors that do no connect to a key that is more recent than current-vio_lag
    if vio:
        for step in range(start_step, end_step+1):
            print(f"step = {step}")
            timestep = timesteps[step]
            max_key = max(timestep.ordering_to_key)

            for i in range(timestep.num_factors - 1, -1, -1):
                remove_flag = False
                factor = timestep.factors[i]

                recent_flag = (max(factor.keys) >= max_key - vio_lag)

                if not recent_flag:
                    remove_flag = True

                if len(factor.keys) > 2 and abs(factor.keys[0] - factor.keys[1] > vio_lag):
                    remove_flag = True

                if remove_flag:
                    timestep.factors[i] = None
                    print(f"Remove {factor} from step {step}")

        # Keys that are older than vio_lag are marginalized
        # Any node that has a connection with a marginalized node is fixed (cannot be relinearized)
        # We should remove factors that span two keys over vio_lag because they will never be processed anyways

        # 1. For each timestep, select relinearized keys based on some strategy
        #    For now, just check delta against relin threshold
        for step in range(start_step, end_step+1):
            timestep = timesteps[step]
            prev_timestep = None if step - 1 < 0 else timesteps[step - 1]

            prev_deltas = [None] if prev_timestep is None else deepcopy(prev_timestep.deltas)

            prev_num_keys = 1 if prev_timestep is None else prev_timestep.num_keys

            # Do not count keys that are older than some age (marginalize)
            # This assumes all factors that span over some time are alrady removed
            
            timestep.relin_keys = []
            min_age = timestep.step - vio_lag
            for key in range(1, prev_num_keys):
                delta_norm = np.linalg.norm(prev_deltas[key], ord=np.inf)
                if key >= min_age and delta_norm > relin_thresh:
                    timestep.relin_keys.append(key)

            timestep.new_keys = []
            for key in range(prev_num_keys, timestep.num_keys):
                timestep.new_keys.append(key)

            print(prev_num_keys)
            print(timestep.relin_keys)
            print(timestep.new_keys)
            print(timestep.is_reconstruct)

        # Print out all the timesteps in header form
        prev_matrix = np.zeros((0, 0))
        for step in range(start_step, end_step+1):
            outfile = f"{outdir}/step-{step}.h"
            timestep = timesteps[step]

            with open(outfile, "w") as fout:
                write_include_guard(fout)

                prev_matrix = timestep.print_header_vio(fout, vio_lag, prev_matrix)

        meta_header = f"{outdir}/vio_dataset.h"
        stepfile_format = "step-{}.h"
        with open(meta_header, "w") as fout:
            write_include_guard(fout)

            Timestep.print_metadata_vio(fout, timesteps, stepfile_format)

    else: # Not vio
        # For each timestep
        # 1. Go through all relin and new keys and mark all the marked cliques and fixed cliques
        #    A clique is marked if it is directly changed. A clique is fixed if its contribution affects a marked clique
        # 2. Determine if a timestep should perform edit or reconstruct
        #    a. If the ordering of keys is different from the previous timestep, it must be a reorder
        #    b. For now, keep everything reconstruct
        # 3. Print out the timestep
        #    a. If reconstructing, print out all the [A B] matrix of each fixed and marked clique

        # This is for full dataset eval
        for step in range(start_step, end_step + 1):
            if step % period != 0:
                print(f"step {step} is set to None")
                timesteps[step] = None

        for step, timestep in enumerate(timesteps):
            if timestep is None:
                continue

            timestep.mark_cliques()
            timestep.fix_cliques()

            prev_timestep = None if step - 1 < 0 else timesteps[step - 1]

            # Determine if edit or reconstruct
            timestep.is_reconstruct = True
            if prev_timestep is None:
                timestep.is_reconstruct = True
            elif not timestep.same_ordering(prev_timestep):
                timestep.is_reconstruct = True  
                timestep.is_reorder = True
            else:
                timestep.is_reconstruct = True
            

            outfile = f"{outdir}/step-{step}.h"
            with open(outfile, "w") as fout:
                write_include_guard(fout)
                if timestep.is_reconstruct:
                    timestep.print_header_reconstruct(fout)


        meta_header = f"{outdir}/incremental_dataset.h"
        stepfile_format = "step-{}.h"
        with open(meta_header, "w") as fout:
            write_include_guard(fout)

            Timestep.print_metadata_incremental(fout, timesteps, stepfile_format)

        # print(f"Write to header {outfile}")




                



