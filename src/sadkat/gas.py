# + ignore="True"
from sadkat.solutes import *
# -

# ## 2.3. Environmental conditions
#
# Data structures for describing the environment (i.e. the gas phase surrounding the droplet):

# +
@dataclass
class Environment:
    """Class to conveniently store all parameters needed to describe the surrounding gas together."""
    solvent: object
    molar_mass: float               # g/mol
    pressure: float                 # Pa
    temperature: float              # K
    relative_humidity: float        # unitless, bounded in [0,1]
    specific_heat_capacity: float   # J/kg/K
    thermal_conductivity: float     # J/s/m/K
    dynamic_viscosity: float        # kg/m/s
    velocity: np.array=np.zeros(3)  # m/s

    @property
    def density(self):
        """Density of gas assuming ideal gas law in kg/m^3."""
        return (1e-3*self.molar_mass) * self.pressure / (gas_constant * self.temperature) # kg/m^3

    @property
    def vapour_pressure(self):
        """Vapour pressure in Pascals."""
        return self.relative_humidity * self.solvent.equilibrium_vapour_pressure(self.temperature)
# -

# Specific parameterisations for the Earth's atmosphere (the main environment people will likely be considering):

# +
def thermal_conductivity_air(T):
    """Thermal conductivity of air in J/s/m/K.

    Source: Stephan and Laesecke, J. Phys. Chem Ref. Data, 14, 1 (1985)
    """
    coefficients = np.array([  33.9729025, -164.702679,  262.108546 , - 21.5346955,
                             -443.455815 ,  607.339582 , -368.790121,  111.296674 ,
                             - 13.4122465])
    exponents = np.linspace(-1, 5/3, len(coefficients)).reshape(-1,1)

    # Fit is for reduced quantities, so convert via:
    Tc = 132.52 # K
    Tr = T / Tc
    Lam = 4.358e-3 # J/s/m/K

    return Lam * coefficients.dot(Tr**exponents)

specific_heat_capacity_air = 1.006e3 # J/kg/K
dynamic_viscosity_air = 1.81e-5 # kg/m/s

def Atmosphere(temperature,
               relative_humidity=0,
               pressure=standard_atmospheric_pressure,
               velocity=np.zeros(3)):
    """Set up conditions for Earth's atmosphere.

    Args:
        temperature: room temperature (K)
        relative_humidity: RH of water vapour (default=0 for dry air).
        pressure: room pressure (default=1 atm) (Pa)
        velocity: velocity in (m/s, dimensional vector)
    Returns:
        Environment object describing room conditions.
    """
    vapour_pressure_water = relative_humidity * Water.equilibrium_vapour_pressure(temperature)
    mole_fraction_water = vapour_pressure_water / pressure
    molar_mass = (1-mole_fraction_water) * molar_mass_dry_air + mole_fraction_water * Water.molar_mass

    return Environment(Water, molar_mass, pressure, temperature, relative_humidity, specific_heat_capacity_air, thermal_conductivity_air(temperature), dynamic_viscosity_air, velocity)
# -

# Sanity check parameterisations of Earth's atmosphere by plotting key quantities below:

if __name__ == '__main__':
    T_C = np.linspace(0, 100, 101)
    fig, (ax1, ax2) = plt.subplots(ncols=2, dpi = resolution)

    first = True
    for RH in np.linspace(0, 1, 6):
        label=('%.1f' % RH)
        if first:
            label = 'RH=%s' % label
            first = False

        rho = [Atmosphere(T+T_freezing, RH).density for T in T_C]
        ax1.plot(T_C, rho, label=label)

    ax1.legend(loc='best')
    ax1.set_xlabel('T (℃)')
    ax1.set_ylabel('density of air (kg/m$^3$)')

    T = np.linspace(150, 1400, 1401)
    ax2.plot(T, 1e3*thermal_conductivity_air(T))
    ax2.set_xlim([0, 1400])
    ax2.set_ylim([0, 100])
    ax2.set_xlabel('T (K)')
    ax2.set_ylabel('thermal conductivity of air (mW/m/K)')

    from matplotlib.collections import LineCollection

    T_C = np.linspace(0, 100, 101)
    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize = figure_size, dpi = resolution)

    cmap = mpl.cm.plasma_r

    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    colors = cmap(np.linspace(0, 1, 101))

    #first = True
    for RH, color in zip(np.linspace(0, 1, 101), colors):
        #label=('%.1f' % RH)
        #if first:
            #label = 'RH=%s' % label
            #first = False

        rho = [Atmosphere(T+T_freezing, RH).density for T in T_C]
        lwidths=(T_C/50)[:-1]
        points = np.array([T_C+T_freezing, rho]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments, linewidths=lwidths,color=color)
        ax1.add_collection(lc)
        ax1.plot(T_C + T_freezing, rho, color = color, lw = 0 , label=label)

    #ax1.legend(loc='best')
    ax1.set_xlabel('T / K')
    ax1.set_ylabel('Density / kgm$^{-3}$')
    fig.colorbar(mpl.cm.ScalarMappable(mpl.colors.Normalize(vmin=0, vmax=100), cmap), ax = ax1,
                 orientation = 'horizontal', label = 'Relative Humidity / %')


    ax2.plot(T_C + T_freezing, 1e3*thermal_conductivity_air(T_C + T_freezing), lw = 3, color = 'k')
    ax2.set_xlim([min(T_C + T_freezing), max(T_C + T_freezing)])
    #ax2.set_ylim([0, 100])
    ax2.set_xlabel('T / K')
    ax2.set_ylabel('Thermal Conductivity / mWm$^{-1}$K$^{-1}$')

    from matplotlib.collections import LineCollection

    T_C = np.linspace(0, 100, 101)
    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize = figure_size, dpi = resolution)

    cmap = mpl.cm.plasma_r

    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    colors = cmap(np.linspace(0, 1, 101))

    #first = True
    for RH, color in zip(np.linspace(0, 1, 101), colors):
        #label=('%.1f' % RH)
        #if first:
            #label = 'RH=%s' % label
            #first = False

        rho = [Atmosphere(T+T_freezing, RH).density for T in T_C]
        lwidths=(T_C/50)[:-1]
        points = np.array([T_C+T_freezing, rho]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments, linewidths=lwidths,color=color)
        ax1.add_collection(lc)
        ax1.plot(T_C + T_freezing, rho, color = color, lw = 0 , label=label)

    #ax1.legend(loc='best')
    ax1.set_xlabel('T$_{gas}$ / K')
    ax1.set_ylabel('Density / kgm$^{-3}$')
    fig.colorbar(mpl.cm.ScalarMappable(mpl.colors.Normalize(vmin=0, vmax=100), cmap), ax = ax1,
                 orientation = 'horizontal', label = 'Relative Humidity / %')


    ax2.plot(T_C + T_freezing, 1e3*thermal_conductivity_air(T_C + T_freezing), lw = 3, color = 'k')
    ax2.set_xlim([min(T_C + T_freezing), max(T_C + T_freezing)])
    #ax2.set_ylim([0, 100])
    ax2.set_xlabel('T$_{gas}$ / K')
    ax2.set_ylabel('Thermal Conductivity / mWm$^{-1}$K$^{-1}$')

    plt.show()
