# # 1. Preamble: standard modules and parameters

# ## 1.1. Importing modules

# +
from dataclasses import dataclass, field
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

# +
from IPython import get_ipython
running_as_jupyter_notebook = 'zmqshell' in str(type(get_ipython()))
if running_as_jupyter_notebook:
    # Style settings to make figures legible in Jupyter notebook
    plt.style.use('default')                # override any previous changes (from e.g. user's matplotlibrc file)
    plt.style.use('notebookStyle.mplstyle') # apply new settings
# -
