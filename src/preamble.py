# # 1. Preamble: standard modules and parameters

# ## 1.1. Importing modules

# +
from dataclasses import dataclass
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

import matplotlib as mpl, matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

import chemicals
import os

# Set up the notebook so that the plots look about the right size.
# %matplotlib inline
#plt.rcParams['figure.figsize'] = [9.5, 6]

figure_width = 360 / 25.4 #conversion to mm is 25.4
figure_height = 222.48 / 25.4 #conversion to mm is 25.4
figure_size = (figure_width, figure_height)

resolution = 600 #dpi
tick_size = 18
fontlabel_size = 18

params = {
    'lines.markersize' : 2,
    'axes.labelsize': fontlabel_size,
    'legend.fontsize': fontlabel_size,
    'xtick.labelsize': tick_size,
    'ytick.labelsize': tick_size,
    'figure.figsize': figure_size,
    'xtick.direction':     'in',     # direction: {in, out, inout}
    'ytick.direction':     'in',     # direction: {in, out, inout}
    'axes.spines.top': False,
    'axes.spines.right': False,
    'xtick.major.pad':  8,
    'ytick.major.pad':  8,
    'axes.linewidth' : 1.2,
    'font.family' : 'serif'}

plt.rcParams.update(params)
# -

# ## 1.2. Standard parameters

# +
T_freezing = 273.15                         # K
standard_temperature = T_freezing + 20      # K
body_temperature = T_freezing + 37          # K
standard_atmospheric_pressure = 101325      # Pa
gas_constant = 8.314472                     # J/K/mol
stefan_boltzmann_constant = 5.670374419e-8  # W/m^2/K^4

molar_mass_dry_air = 28.9647                # g/mol

gravitational_acceleration = 9.80665        # m/s^2
# -

