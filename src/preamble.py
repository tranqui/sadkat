# # 1. Preamble: standard modules and parameters

# ## 1.1. Importing modules

# +
from dataclasses import dataclass
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Set up the notebook so that the plots look about the right size.
# %matplotlib inline
plt.rcParams['figure.figsize'] = [9.5, 6]
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

