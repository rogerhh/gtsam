import yaml
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.stats import special_ortho_group
import sys
from optparse import OptionParser
import math
from copy import deepcopy
from utils import *
import glob
import os
from subprocess import Popen, PIPE, DEVNULL, CalledProcessError
from matplotlib.animation import FuncAnimation
import threading

def get_yes_no_input(prompt):
    while True:
        user_input = input(prompt).strip().lower()
        if user_input in {'y', 'n'}:
            return user_input == 'y'
        print("Invalid input. Please enter 'Y' or 'N'.")

def run_system_cmd(cmd):
    print(cmd)
    with Popen(cmd, shell=True, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
        for line in p.stdout:
            print(line, end='') # process line here

        p.wait()

        if p.returncode != 0:
            print(p.stderr)
            raise CalledProcessError(p.returncode, p.args)

def convert_to_numeric(variable):
    if isinstance(variable, (int, float)):
        return variable
    elif variable[-1] == 'M':
        return int(variable[:-1]) * (2 ** 20)
    elif variable[-1] == 'K':
        return int(variable[:-1]) * (2 ** 10)
    else:
        raise ValueError("Variable format not recognized")

def run_config(cmd, pipe, traj, chis, lock, latencies=None, alg=None):
    print("In run_config")
    print(cmd, pipe)
    if not os.path.exists(pipe):
        os.mkfifo(pipe)

    p = Popen(cmd, shell=True, stdout=DEVNULL, bufsize=1, universal_newlines=True)

    last_received_chi = 0

    with open(pipe, "r") as fifo:
        # for line in p.stdout:
        #     print(line, end='') # process line here
        step = 0
        while True:
            # print("reading")
            # for line in p.stdout:
            #     print(line, end='')
            line = fifo.readline()

            if line == "":
                print("end of pipe")
                break

            if "start" in line:
                step = 0
            elif "end" in line:
                line = fifo.readline()
                assert("chi2" in line)

                chi2 = float(line.split()[2])
                chi2_copy = chi2

                if last_received_chi == 0:
                    chi2 = 0

                last_received_chi = chi2_copy

                chis.append(chi2)

                if alg == "ra":
                    for _ in range(10):
                        line = fifo.readline()
                        if "totalCost:" in line:
                            latency = float(line.split()[-1])
                            latencies.append(latency)
                            break
                else:
                    latencies.append(0)

            else:
                arr = line.split()
                M = np.zeros((3, 4));
                for i in range(3):
                    for j in range(4):
                        M[i, j] = float(arr[i * 4 + j])
                x, y, z = M[0, 3], M[1, 3], M[2, 3]

                lock.acquire()
                if len(traj["x"]) <= step:
                    traj["x"].append(x)
                    traj["y"].append(y)
                    traj["z"].append(z)
                else:
                    traj["x"][step] = x
                    traj["y"][step] = y
                    traj["z"][step] = z

                lock.release()

                step += 1

    os.remove(pipe)

class Updater:
    def __init__(self, labels, axes, lines, trajs, chis, locks, chi_axes=None, chi_lines=None, latencies=None, latency_axes=None, latency_lines=None):
        self.labels = labels
        self.axes = axes
        self.trajs = trajs
        self.locks = locks
        self.lines = lines
        self.chis = chis
        self.chi_axes = chi_axes
        self.chi_lines = chi_lines
        self.latency = latencies
        self.latency_axes = latency_axes
        self.latency_lines = latency_lines

    def update(self, frame):
        for i in range(len(self.locks)):
            lock = self.locks[i]
            ax = self.axes[i]
            traj = self.trajs[i]
            line = self.lines[i]
            chi = self.chis[i]
            chi_ax = self.chi_axes[i]
            chi_line = self.chi_lines[i]

            latency = self.latency[i]
            latency_ax = self.latency_axes[i]
            latency_line = self.latency_lines[i]

            if lock.acquire(False):

                x = traj["x"]
                y = traj["y"]
                z = traj["z"]
                
                line.set_xdata(x)
                line.set_ydata(y)
                line.set_3d_properties(z)

                chi_line.set_xdata(range(len(chi)))
                chi_line.set_ydata(chi)

                latency_line.set_xdata(range(len(latency)))
                latency_line.set_ydata(latency)

                lock.release()

                ax.set_xlim(min(x) - 5, max(x) + 5)
                ax.set_ylim(min(y) - 5, max(y) + 5)
                ax.set_zlim(min(z) - 5, max(z) + 5)

                chi_ax.set_xlim(-1, len(chi) + 5)
                chi_ax.set_ylim(1e-10, max(chi) * 5)

                chi_ax.legend([labels[i]])

                latency_ax.set_xlim(-1, len(latency) + 5)
                latency_ax.set_ylim(0, max(latency) + 5)




if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("--config", dest="config", 
                      default="config_demo.yaml", help="The path to the config file")
    (options, args) = parser.parse_args()

    with open(options.config, "r") as config_fin:
        config = yaml.safe_load(config_fin)

    #######################################
    # Metadata
    #######################################
    scriptdir = os.path.abspath(config["metadata"]["scriptdir"])
    headerdir = os.path.abspath(config["metadata"]["headerdir"])
    builddir = os.path.abspath(config["metadata"]["builddir"])
    recompile = config["metadata"]["recompile"]
    make_flags = config["metadata"]["make_flags"]

    binary_str = ""
    for binary in config["metadata"]["target_binary"]:
        binary_str += f"{binary} "

    if recompile:
        run_system_cmd(f"cd {builddir} && cmake .. && make {make_flags} {binary_str}")

    dataset = config["target_dataset"]

    labels = []
    cmds = []
    pipes = []
    algs = []
    latency_targets = [3.33 * 10^7 for _ in cmds]

    for config_name in config["target_configs"]:
        labels.append(config_name)
        algs.append(config[config_name]["alg"])

        d = config[config_name]
        d.update(config["dataset"][dataset])

        is3D = d["is3D"]
        dataset_path = d["path"]
        relin_thresh = d["relin_thresh"]
        start_step = d["start_step"]
        end_step = d["end_step"] + 1
        period = d["period"]
        pipe_name = f"/tmp/{config_name}_pipe"
        print(d)

        if d["alg"] == "ra":
            num_threads = d["num_threads"]
            latency = d["latency"]

            exe = f"{builddir}/timing/" + ("testGtsamIncremental3D-ra" if is3D else "testGtsamIncremental-ra")

            print(f"Running config {config_name}")
            cmd = f"{exe} -f {dataset_path} \
                          --num_steps {end_step} \
                          --relin_thresh {relin_thresh} \
                          --num_threads {num_threads} \
                          {'--noise_format {}'.format(noise_format) if not is3D else ''} \
                          --print_frequency 1 \
                          --demo_pipe {pipe_name} \
                          "
            # run_system_cmd(cmd)

        if d["alg"] == "vio":
            max_iter = d["max_iter"]
            vio_lag = d["vio_lag"]
            run_lc = d["run_lc"]
            lc_period = d["lc_period"]

            exe = f"{builddir}/timing/" + ("testGtsamIncremental3D-vio" if is3D else "testGtsamIncremental-vio")

            cmd = f"{exe} -f {dataset_path} \
                          --num_steps {end_step} \
                          --relin_thresh {relin_thresh} \
                          -m {max_iter} \
                          -e 0 \
                          -d 0 \
                          --vio_lag {vio_lag} \
                          {'--run_lc' if run_lc else ''} \
                          {'--lc_period {}'.format(lc_period) if run_lc else ''} \
                          {'--noise_format {}'.format(noise_format) if not is3D else ''} \
                          --print_frequency 1 \
                          --demo_pipe {pipe_name} \
                          "
        cmds.append(cmd)
        pipes.append(pipe_name)


    threads = []
    trajs = [{"x":[], "y":[], "z":[]} for _ in cmds]
    chis = [[] for _ in cmds]
    latencies = [[] for _ in cmds]
    locks = [threading.Lock() for _ in cmds]

    fig = plt.figure(figsize=(20, 30))

    num_plots = len(cmds)
    axes = [fig.add_subplot(3, num_plots, i+1, projection="3d") for i in range(num_plots)]
    lines = [ax.plot([], [], [])[0] for ax in axes]

    chi_axes = [fig.add_subplot(3, num_plots, num_plots + i+1) for i in range(num_plots)]
    chi_lines = [chi_ax.plot([], [])[0] for chi_ax in chi_axes]

    latency_axes = [fig.add_subplot(3, num_plots, 2 * num_plots + i+1) for i in range(num_plots)]
    latency_lines = [latency_ax.plot([], [])[0] for latency_ax in latency_axes]

    for chi_ax in chi_axes:
        chi_ax.set_yscale("log")

    updater = Updater(labels, axes, lines, trajs, chis, locks, chi_axes=chi_axes, chi_lines=chi_lines, latencies=latencies, latency_axes=latency_axes, latency_lines=latency_lines)

    for cmd, pipe, traj, chi, lock, alg, latency in zip(cmds, pipes, trajs, chis, locks, algs, latencies):
        
        thread = threading.Thread(target=run_config, args=(cmd, pipe, traj, chi, lock, latency, alg))
        thread.start()
        threads.append(thread)

    animation = FuncAnimation(fig, updater.update, interval=200)

    plt.show()

    for thread in threads:
        thread.join()



