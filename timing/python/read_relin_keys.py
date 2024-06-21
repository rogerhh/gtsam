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
    parser.add_option("--outfile", dest="outfile", 
                      default=None, help="The path to the output file")
    parser.add_option("--start_step", dest="start_step", type="int",
                      default=2, help="Starting step, inclusive")
    parser.add_option("--end_step", dest="end_step", type="int",
                      default=30, help="Ending step, inclusive")
    (options, args) = parser.parse_args()

    indir = options.indir
    outfile = options.outfile
    start_step = options.start_step
    end_step = options.end_step
    period = 1

    files = glob.glob(os.path.join(indir, "step-*.out"))

    timesteps = [None for _ in range(0, end_step+1)]

    infiles = [None for _ in range(0, end_step+1)]

    infile_count = 0

    for step in range(start_step, end_step+1):
        if step % period != 0:
            continue

        infile_pattern = f"step-{step}.out"
        infile = None

        for filename in files:
            if infile_pattern in filename:
                infile = filename
                infiles[step] = infile
                infile_count += 1
                break

    
    with open(outfile, "w") as fout:
        fout.write(f"{infile_count}\n")
        for step in range(start_step, end_step+1):
            if step % period != 0:
                continue

            infile = infiles[step]

            if infile is None:
                print(f"{step} file not found!", step)
                assert(0)

            with open(infile, "r") as fin:
                read_until(fin, "relin keys")

                line = fin.readline()

                num_relin_keys = int(line)

                relin_keys = []

                fout.write(f"{step} {num_relin_keys}\n")

                for i in range(num_relin_keys):
                    relin_key = int(fin.readline())
                    fout.write(f"{relin_key} ")
                fout.write("\n")

