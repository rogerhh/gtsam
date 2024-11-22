import sys
from optparse import OptionParser
from copy import deepcopy

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("--infile", dest="infile",
                      default=None, help="input file path")
    parser.add_option("--outfile", dest="outfile",
                      default=None, help="output file path")
    (options, args) = parser.parse_args()

    outlines = []
    with open(options.infile, "r") as fin:
        lines = fin.readlines()
        for line in lines:
            if "EDGE" not in line:
                outlines.append(line)
            else:
                arr = line.split()
                key1 = int(arr[1])
                key2 = int(arr[2])

                if abs(key1 - key2) >= 4 and abs(key1 - key2) <= 8:
                    continue
                else:
                    outlines.append(line)

    with open(options.outfile, "w") as fout:
        for line in outlines:
            fout.write(line)



