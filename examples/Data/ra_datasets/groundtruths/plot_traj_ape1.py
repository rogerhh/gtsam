import sys
from optparse import OptionParser
from copy import deepcopy
import matplotlib.pyplot as plt

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("--infiles", dest="infiles",
                      default=None, help="input file paths")
    parser.add_option("--gt_infile", dest="gt_infile",
                      default=None, help="ground truth file path")
    parser.add_option("--labels", dest="labels",
                      default=None, help="labels")
    parser.add_option("--start_step", dest="start_step", type="int",
                      default=2, help="input file path")
    parser.add_option("--outfile", dest="outfile",
                      default=None, help="output file")
    (options, args) = parser.parse_args()

    infiles = options.infiles.split()
    gt_infile = options.gt_infile
    labels = options.labels.split()
    start_step = options.start_step

    plt.rcParams.update({'font.size': 42, 'axes.titlesize': "large", 'axes.labelsize': 'large'})
    fontsize = 42

    t = []
    rmse = []
    max_error = []
    min_error = []

    figure = plt.figure(figsize=(30, 9))

    ax1 = plt.subplot(1, 2, 2)
    ax2 = plt.subplot(1, 2, 1)

    gt_max = []
    gt_rmse = []

    # with open(gt_infile, "r") as fin:
    #     while True:
    #         line = fin.readline()
    #         if not line:
    #             break

    #         arr = line.split()

    #         if "not aligned" in line:
    #             # t[-1].append(start_step)
    #             # start_step += 1
    #             pass

    #         if "max" in line:
    #             e = float(arr[1])
    #             gt_max.append(e)
    #         if "rmse" in line:
    #             e = float(arr[1])
    #             gt_rmse.append(e)


    for i, infile in enumerate(infiles):
        t.append([])
        rmse.append([])
        max_error.append([])
        start_step = options.start_step
        with open(infile, "r") as fin:
            index = -1
            while True:
                line = fin.readline()
                if not line:
                    break

                arr = line.split()

                if "not aligned" in line:
                    t[-1].append(start_step)
                    start_step += 1
                    index += 1

                if "max" in line:
                    e = float(arr[1])
                    max_error[-1].append(e)
                if "rmse" in line:
                    e = float(arr[1])
                    rmse[-1].append(e)

        ax1.plot(t[-1], rmse[-1], label=f"{labels[i]}", linewidth=5)
        ax2.plot(t[-1], max_error[-1], label=f"{labels[i]}", linewidth=5)
    ax1.set_xlabel("Step", fontsize=fontsize)
    ax1.set_ylabel("Error (m)", fontsize=fontsize)
    ax1.set_title("RMSE", fontsize=fontsize)

    tick_size = 42
    ax1.tick_params(axis='both', which='major', labelsize=tick_size)
    # ax1.set_ylim([-5, 15])
    # ax1.set_yscale("log")

    ax2.set_title("MAX", fontsize=fontsize)
    ax2.set_xlabel("Step")
    ax2.set_ylabel("Error (m)", fontsize=fontsize)
    # ax2.set_yscale("log")
    ax2.tick_params(axis='both', which='major', labelsize=tick_size)
    # ax2.set_ylim([-5, 15])

    figure.tight_layout()


    legend_fontsize = 30
    ax1.legend(loc="upper left", fontsize=legend_fontsize)

    plt.savefig(options.outfile, bbox_inches='tight')



