import numpy as np
import scipy
from scipy.stats import special_ortho_group
import sys
from optparse import OptionParser
import math
from copy import deepcopy

if __name__ == "__main__":
    parser = OptionParser()

    parser.add_option("--delta_dir", dest="indir", type=str,
                      help="Directory containing delta files")
    parser.add_option("--output_dir", dest="outdir", type=str,
                      help="Directory to write output header files")
    parser.add_option("--dim", dest="dim", type=int,
                      help="Dimension of the variables")

    (options, args) = parser.parse_args()

    prev_deltas = []

    step = 2
    while True:
        filename = f"{options.indir}/step-{step}_delta.out"

        # Check if file exists
        try:
            with open(filename, "r") as fin:
                lines = fin.readlines()
                n = int(lines[1].strip())
                cur_deltas = []
                indices = []
                for i in range(n):
                    delta = np.array([float(x) for x in lines[i + 2].split()])
                    dim = len(delta)
                    cur_deltas.append(delta)
                    if len(prev_deltas) > i and np.allclose(delta, prev_deltas[i], rtol=1e-4):
                        pass
                    else:
                        indices.append(i)

                l = len(indices)

                # Generate header file
                header = f"{options.outdir}/step-{step}_delta.h"
                with open(header, "w") as fout:
                    fout.write(f"#ifndef STEP_{step}_DELTA_H\n")
                    fout.write(f"#define STEP_{step}_DELTA_H\n\n")
                    fout.write(f"const int INJECT_DELTA_step_{step}_delta_n = {n};\n")
                    fout.write(f"const int INJECT_DELTA_step_{step}_indices_n = {l};\n")
                    fout.write(f"const int INJECT_DELTA_step_{step}_indices[{l}] = {{")
                    for i in range(len(indices)):
                        fout.write(f"{indices[i]}")
                        if i == len(indices) - 1:
                            fout.write("};\n\n")
                        else:
                            fout.write(", ")
                    fout.write(f"const double INJECT_DELTA_step_{step}_delta_data[{dim * l}] = {{\n")
                    for i in range(l):
                        delta = cur_deltas[indices[i]]
                        for j in range(dim):
                            fout.write(f"{delta[j]}")
                            if i == l - 1 and j == dim - 1:
                                fout.write("};\n\n")
                            else:
                                fout.write(", ")
                        fout.write("\n")

                    fout.write("#endif\n")
            step += 1
        except FileNotFoundError:
            step -= 1
            break

    # Generate meta header file
    max_steps = step
    meta_header = f"{options.outdir}/incremental_delta.h"
    with open(meta_header, "w") as fout:
        fout.write("#ifndef INCREMENTAL_DELTA_H\n")
        fout.write("#define INCREMENTAL_DELTA_H\n\n")

        fout.write("#define INJECT_DELTA_HEADER\n\n")

        for i in range(2, max_steps + 1):
            fout.write(f"#include \"step-{i}_delta.h\"\n")
        fout.write("\n")

        fout.write(f"const int INJECT_DELTA_dim = {dim};\n")

        fout.write(f"const int INJECT_DELTA_incremental_delta_steps = {step};\n\n")

        fout.write(f"const int INJECT_DELTA_step_delta_n[{max_steps+1}] = {{\n")
        for i in range(2):
            fout.write(f"0, ")
        for i in range(2, max_steps + 1):
            fout.write(f"INJECT_DELTA_step_{i}_delta_n, ")
        fout.write("};\n\n")

        fout.write(f"const int INJECT_DELTA_step_indices_n[{max_steps+1}] = {{\n")
        for i in range(2):
            fout.write(f"0, ")
        for i in range(2, max_steps + 1):
            fout.write(f"INJECT_DELTA_step_{i}_indices_n, ")
        fout.write("};\n\n")

        fout.write(f"const int* INJECT_DELTA_step_indices[{max_steps+1}] = {{\n")
        for i in range(2):
            fout.write(f"nullptr, ")
        for i in range(2, max_steps + 1):
            fout.write(f"INJECT_DELTA_step_{i}_indices, ")
        fout.write("};\n\n")

        fout.write(f"const double* INJECT_DELTA_step_delta_data[{max_steps+1}] = {{\n")
        for i in range(2):
            fout.write(f"nullptr, ")
        for i in range(2, max_steps + 1):
            fout.write(f"INJECT_DELTA_step_{i}_delta_data, ")
        fout.write("};\n\n")

        fout.write("\n#endif\n")
