"""
Author: Roger Hsiao
main.py: Driver program for the python version of gtsam for testing

"""

import sys
from optparse import OptionParser

if __name__ == "__main__":
    # Parse options
    parser = OptionParser()
    parser.add_option("-f", "--dataset", dest="dataset",
                      default="", help="Dataset file in g2o format.")

    # 1. Read in constraints

    # 2. Initialize X_new

    # 3. If reordering, determine best ordering

    # 4. Relinearize

    # 5. Solve min_x \| Ax - b \|^2_2 in whatever way

    # 6. Compute prediction error, if too high, return to 4.
