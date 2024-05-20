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

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from threading import Thread

# Initialize figure and axis
fig = plt.figure()

# Create empty line objects for each subplot
ax1 = fig.add_subplot(121, projection='3d')
line1, = ax1.plot([], [], [], color='blue')
ax2 = fig.add_subplot(122, projection='3d')
line2, = ax2.plot([], [], [], color='red')

# Function to update plot data for subplot 1
def update_subplot1(frame):
    print(frame)
    x = np.linspace(0, 2 * np.pi, 100)
    y = np.sin(x + frame * 0.1)
    z = np.cos(x + frame * 0.1)
    line1.set_xdata(x)
    line1.set_ydata(y)
    line1.set_3d_properties(z)

    ax1.set_xlim(min(x) - 5, max(x) + 5)
    ax1.set_ylim(min(y) - 5, max(y) + 5)
    ax1.set_zlim(min(z) - 5, max(z) + 5)

    x = np.linspace(0, 2 * np.pi, 100)
    y = np.cos(x + frame * 0.1)
    z = np.sin(x + frame * 0.1)
    line2.set_xdata(x)
    line2.set_ydata(y)
    line2.set_3d_properties(z)

    ax2.set_xlim(min(x) - 5, max(x) + 5)
    ax2.set_ylim(min(y) - 5, max(y) + 5)
    ax2.set_zlim(min(z) - 5, max(z) + 5)

# Function to update plot data for subplot 2
def update_subplot2(frame):
    print(frame)
    x = np.linspace(0, 2 * np.pi, 100)
    y = np.cos(x + frame * 0.1)
    z = np.sin(x + frame * 0.1)
    line2.set_xdata(x)
    line2.set_ydata(y)
    line2.set_3d_properties(z)

ani1 = None
def animate(fig, func):
    global ani1
    ani1 = FuncAnimation(fig, func, frames=range(100), interval=200)

animation = FuncAnimation(fig, update_subplot1, interval=200)

# Show the plot
plt.show()

