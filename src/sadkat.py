# # 1. Parameterising the model

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

# ## 1.3. Solvents

# ### 1.3.1. Data structures for general solvents.

# +
@dataclass
class Solvent:
    """Class to conveniently store all parameters needed to describe a solvent together."""

    molar_mass: float                            # g/mol
    density: object                              # kg/m^3
    specific_heat_capacity: object               # J/kg/K
    specific_latent_heat_vaporisation: object    # J/kg
    equilibrium_vapour_pressure: object          # Pa
    vapour_binary_diffusion_coefficient: object  # m^2/s

@dataclass
class VapourBinaryDiffusionCoefficient:
    """Standard fitting function to describe the temperature dependence of the binary diffusion coeffient
    for evaporated solvent in air.

    The funtional form of this is taken from Xie et al., Indoor Air (2007).
    """
    D_ref: float   # m^2/s
    T_ref: float   # K
    lam: float     # unitless, scaling parameter

    def __call__(self, T):
        """Fit itself.

        Args:
            T: temperature in Kelvin.
        Returns:
            The diffusion coefficient in m^2/s.
        """
        return self.D_ref * (T / self.T_ref)**self.lam


# -

# ### 1.3.2. Properties of pure water

# +
molar_mass_water = 18.01528 # g/mol

def density_water(temperature):
    """Fit for the density of pure water used in J Walker model.
    
    Originally from:

        Wagner and Pruss Addendum to J. Phys. Chem. Ref. Data 16, 893 (1987),
        J. Phys. Chem. Ref. Data, 1993, 22, 783–787)

    Args:
        temperature: temperature in Kelvin.
    Returns:
        The water density in kg/m^3.
    """
    ref_T = 647.096 # K
    ref_density = 322 # kg/m^3
    b1, b2, b3, b4, b5, b6 = 1.99274064, 1.09965342, -0.510839303, -1.75493479, -45.5170352, -674694.45

    theta = temperature / ref_T

    tau = 1 - theta

    density = ref_density * (1 + b1*pow(tau,1/3) + b2*pow(tau,2/3) + b3*pow(tau,5/3) + b4*pow(tau,16/3) + b5*pow(tau,45/3) + b6*pow(tau,110/3))

    return density

# # ? not sure where this is from
specific_heat_capacity_water = lambda T: 4180 # J/kg/K

# Su, PCCP (2018)
specific_latent_heat_water = lambda T: 3.14566e6 - 2361.64 * T # J/kg

def equilibrium_vapour_pressure_water(T):
    """Using the Buck equation (from integrating the Clausius-Clapeyron equation).

    Args:
        T: temperature in Kelvin.
    """
    T_C = T - T_freezing # Celsius
    return 1e3*0.61161 * np.exp((18.678 - (T_C / 234.5)) * (T_C / (257.14 + T_C))) # Pa

Water = Solvent(molar_mass_water,
                density_water,
                specific_heat_capacity_water,
                specific_latent_heat_water,
                equilibrium_vapour_pressure_water,
                VapourBinaryDiffusionCoefficient(0.2190e-4, T_freezing, 1.81))
# -

# Sanity check the water properties by plotting them below:

# +
T_C = np.linspace(0, 100, 101)
fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(ncols=2, nrows=2)

ax1.plot(T_C, Water.density(T_C + T_freezing))
ax1.set_xlabel('T (℃)')
ax1.set_ylabel('density (kg/m3)')

ax2.plot(T_C, 1e-3 * Water.specific_latent_heat_vaporisation(T_C + T_freezing))
ax2.axhline(y=2264.705, ls='dashed')
ax2.set_xlabel('T (℃)')
ax2.set_ylabel('L (kJ/kg)')

ax3.plot(T_C, Water.equilibrium_vapour_pressure(T_C+T_freezing))
ax3.set_xlabel('T (℃)')
ax3.set_ylabel('equilibrium vapour pressure (Pa)')

ax4.plot(T_C, Water.vapour_binary_diffusion_coefficient(T_C+T_freezing))
ax4.set_xlabel('T (℃)')
ax4.set_ylabel('$D_\infty$ (m$^2$/s)')

plt.show()


# -

# ## 1.4. Solutes
#
# Parameterising solvent activity in terms of mass fraction of solute, and its inverse. We define the ideal case (i.e. applying Raoult's law for a mixture) and two methods of parameterising the fits in the non-ideal case:
#
# * Parameterising mfs -> activity (with activity -> mfs obtained by inverting this parameterisation).
# * Parameterising activity -> mfs (with mfs -> activity obtained by inverting this parameterisation).
#
# We provide both options for parameterising the non-ideal case to increase flexibility in defining new solutions.

# +
class IdealSolventActivity:
    """In ideal mixtures the solvent activity obeys Raoult's law."""

    def __init__(self, solution):
        """Construct the activity rule for this mixture.

        Args:
            solution: parameters describing the solution. We need this to calculate the mole fractions
            in the mixture to apply Raoult's law."""
        self.solution = solution

    def solvent_activity_from_mass_fraction_solute(self, mfs):
        """Solvent activity is simply the mole fraction of solvent under Raoult's law.

        Args:
            mfs: the mass fraction of the solute, a number bounded between 0 and 1
        """
        return self.solution.mole_fraction_solvent(mfs)

    def mass_fraction_solute_from_solvent_activity(self, aw):
        """Invert Raoult's law to obtain the mass fraction at a particular solvent activity.

        This is useful for e.g. obtaining the equilibrium mass fraction when the water activity matches
        the ambient conditions.

        Args:
            aw: the solvent activity, a number bounded between 0 and 1.
        """
        x_solvent = aw
        x_solute = 1 - aw
        w_solvent = x_solvent * self.solution.solvent.molar_mass
        w_solute = x_solute * self.solution.molar_mass_solute
        return w_solute / (w_solute + w_solvent)

def invert_fit(y, yfit):
    """Invert a fit y(x) to yield x(y).

    Args:
        y: y value (a scalar).
        yfit: parameters of polynomial fit y(x).
    Returns:
        The value x(y).
    """
    x = (yfit - y).roots
    x = np.real(np.min(x[np.isreal(x)]))
    return x

class ActivityVsMfsParameterisation:
    """Fit of solvent activity for the forward transformation aw = aw(mfs)."""

    def __init__(self, coefficients):
        """Set up the fit for the Taylor expansion

        aw = mfs + b*mfs + c*mfs**2 + ...

        And its inverse.

        Args:
            coefficients: sequence of coefficients a, b, c, ... etc.
        """
        self.solvent_activity_from_mass_fraction_solute = np.poly1d(np.flipud(coefficients))
        self.mass_fraction_solute_from_solvent_activity = np.vectorize(lambda mfs: invert_fit(mfs, self.solvent_activity_from_mass_fraction_solute))

class MfsVsActivityParameterisation:
    """Fit of solvent activity for the backward transformation mfs = mfs(aw)."""

    def __init__(self, coefficients):
        """Set up the fit for the Taylor expansion

        mfs = a + b*aw + c*aw**2 + ...

        And its inverse.

        Args:
            coefficients: sequence of coefficients a, b, c, ... etc.
        """
        self.mass_fraction_solute_from_solvent_activity = np.poly1d(np.flipud(coefficients))
        self.solvent_activity_from_mass_fraction_solute = np.vectorize(lambda aw: invert_fit(aw, self.mass_fraction_solute_from_solvent_activity))


# -

# Check the fit inversions in each kind of parameterisations are working correctly. That is, if we have a parameterisation from mfs -> activity, then this should be consistent with the inverse mapping activity -> mfs. In each plot below the regular and inverted fits should overlap:

# +
# Dummy parameters for the fits.
forward_fit = ActivityVsMfsParameterisation([1, -1, 1, -1])
backward_fit = MfsVsActivityParameterisation([1, -1, 1, -1])

mfs = np.linspace(0, 1, 100)
aw = np.linspace(0, 1, 100)

fig, (ax1, ax2) = plt.subplots(ncols=2)

ax1.plot(mfs, forward_fit.solvent_activity_from_mass_fraction_solute(mfs), label='regular fit')
ax1.plot(forward_fit.mass_fraction_solute_from_solvent_activity(aw), aw, label='inverted fit')
ax1.legend(loc='best')
ax1.set_title('forward parameterisation')

ax2.plot(backward_fit.mass_fraction_solute_from_solvent_activity(aw), aw, label='regular fit')
ax2.plot(mfs, backward_fit.solvent_activity_from_mass_fraction_solute(mfs), label='inverted fit')
ax2.legend(loc='best')
ax2.set_title('backward parameterisation')
plt.show()


# -

# Data structures for parameterising solutes:

# +
@dataclass
class Solution:
    """Class to conveniently store all parameters needed to describe a solute together."""

    solvent: object
    molar_mass_solute: float           # g/mol
    num_ions: int                      # unitless
    solubility_limit: float            # kg per kg solvent
    solid_density: float               # kg/m^3
    density: object                    # kg/m^3
    solvent_activity_fit: object=None  # a fitting function for solvent activity (unitless) in terms of 
                                       # the mass fraction of solute, or None (the default) to assume
                                       # Raoult's law for ideal mixtures

    def mole_fraction_solute(self, mass_fraction_solute):
        """Mole fraction from mass fraction."""
        w_solute = mass_fraction_solute / self.molar_mass_solute
        w_solvent = (1 - mass_fraction_solute) / self.solvent.molar_mass
        return w_solute / (w_solute + w_solvent)

    def mole_fraction_solvent(self, mass_fraction_solute):
        """Mole fraction of *solvent* from mass fraction of *solute*."""
        return 1 - self.mole_fraction_solute(mass_fraction_solute)

    def solvent_activity(self, mass_fraction_solute):
        """Solvent activity describes deviation in vapour pressure from pure solvent.

        If no solvent activity fit has been specified, we will default to the activity
        of ideal mixtures (Raoult's law).

        Args:
            mass_fraction_solute: as described (0=pure solvent, 1=pure solute)
        """

        if self.solvent_activity_fit is None:
            self.solvent_activity_fit = IdealSolventActivity(self)

        return self.solvent_activity_fit.solvent_activity_from_mass_fraction_solute(mass_fraction_solute)

    def mass_fraction_solute_from_solvent_activity(self, solvent_activity):
        """Invert solvent activity to find the mass fraction of the solute at a particular solvent activity.

        This is useful for determining the equilibrium state of a droplet by matching vapour pressures
        at the liquid-gas boundary.

        Args:
            solvent_activity: the activity to invert at (should be bounded between 0 and 1).
        """

        if self.solvent_activity_fit is None:
            self.solvent_activity_fit = IdealSolventActivity(self)

        return self.solvent_activity_fit.mass_fraction_solute_from_solvent_activity(solvent_activity)

class DensityVsMassFractionFit:
    """Fitting function preferred by Bristol Aerosol Research Centre.

    This function expands the density in terms of half-integer powers of density, i.e.

        density = a + b*mfs^0.5 + c*mfs + d*mfs^1.5 + ...

    are the first few terms.
    """

    def __init__(self, coefficients):
        """Construct the fitting function from a set of coefficients.
        
        Args:
            coefficients: the set of coefficients in the equation above (i.e. the sequence [a,b,c,d,...] etc)
        """
        self.density_vs_sqrt_mfs_fit = np.poly1d(np.flipud(coefficients))

    def __call__(self, mfs):
        """The fitting function itself.

        Args:
            mfs: mass fraction of solute
        Returns:
            The density at this value of mfs.
        """
        return self.density_vs_sqrt_mfs_fit(np.sqrt(mfs))

class VolumeAdditivityFit:
    """Density fit for ideal mixtures assumes the volumes of the solvent/solute combine additively.

    This form will be less accurate than empirical fits, but it is convenient for some theoretical
    calculations.
    """

    def __init__(self, pure_solvent_density, pure_solute_density):
        """The fit is parameterised by the densities of the pure constituents.

        Args:
            pure_solvent_density: density of the pure solvent component (kg/m^3).
            pure_solute_density: density of the pure solute component (kg/m^3).
        """
        self.pure_solvent_density = pure_solvent_density
        self.pure_solute_density = pure_solute_density

        self.mass_difference_parameter = \
            (pure_solvent_density - pure_solute_density) / \
            (pure_solvent_density * pure_solute_density)

    def __call__(self, mfs):
        """The fitting function itself.

        Args:
            mfs: mass fraction of solute
        Returns:
            The density at this value of mfs.
        """
        return 1 / (mfs / self.pure_solute_density + (1-mfs) / self.pure_solvent_density)


# -

#  Parameterisations of some common solutes dissolved in water:

# +
aqueous_NaCl = Solution(Water, 58.44, 2, 0.3, 2170,
                        DensityVsMassFractionFit([998.2 , -55.33776, 1326.69542, -2131.05669, 2895.88613, -940.62808]))

# Warning: besides the density fits, the parameters below for DLF/AS are just dummy placeholder values
#          as far as I am aware.

aqueous_DLF  = Solution(Water,  1.00, 2, 1.0, 1700,
                        DensityVsMassFractionFit([995.78, 262.92   , -606.15   ,  1135.53   ,    0      ,    0]))

aqueous_AS   = Solution(Water,  1.00, 2, 1.0, 1200,
                        DensityVsMassFractionFit([997   , -43.88148,  397.75347,  -100.99474,    0      ,    0]))

all_solutions = {'NaCl in water': aqueous_NaCl,
                 'DLF in water': aqueous_DLF,
                 'AS in water': aqueous_AS}

# Identical versions of the above solutions but with a simpler density parameterisation (volume additivity):
from copy import copy
all_solutions_volume_additivity = {}
for label, solution in all_solutions.items():
    additive_solution = copy(solution)
    rho1, rho2 = solution.density(0), solution.solid_density
    additive_solution.density = VolumeAdditivityFit(rho1, rho2)
    all_solutions_volume_additivity[label] = additive_solution
# -

# Sanity check the parameterisations of the solutions by plotting some of their properties below:

# +
mfs = np.linspace(0, 1, 100)
fig, (ax1, ax2) = plt.subplots(ncols=2)

for label, solute in all_solutions.items():
    additive_solute = all_solutions_volume_additivity[label]

    pl, = ax1.plot(mfs, solute.density(mfs), label=label)
    ax1.plot(mfs, additive_solute.density(mfs), '--', c=pl.get_color(), label=('%s (volume additivity)' % label))
    ax1.plot(1, solute.solid_density, 'o', mfc='None', c=pl.get_color())

ax1.set_xlabel('MFS')
ax1.set_ylabel('density (kg/m$^3$)')
ax1.legend(loc='best')

aw = np.linspace(0, 1, 100)
for label, solute in all_solutions.items():
    pl, = ax2.plot(mfs, solute.solvent_activity(mfs), label=label)
    ax2.plot(solute.mass_fraction_solute_from_solvent_activity(aw), aw, '--', c=pl.get_color(), label=('%s (inverted fit)' % label))

ax2.set_xlabel('MFS')
ax2.set_ylabel('solvent activity')
ax2.legend(loc='best')

plt.show()


# -

# ## 1.5. Environmental Conditions
#
# Data structures for describing the environment (i.e. the gas phase surrounding the droplet):

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

# +
T_C = np.linspace(0, 100, 101)
fig, (ax1, ax2) = plt.subplots(ncols=2)

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

plt.show()


# -

# # 2. Defining the droplet and its evolution
#
# NB: See the appendix at the bottom of this notebook for the complete set of equations used to describe the droplet.

@dataclass
class UniformDroplet:
    """This class completely describes the state of the droplet during its evolution.

    The solute is assumed to be distributed uniformly throughout the droplet.
    """

    solution: object
    environment: object
    mass_solute: float              # kg
    mass_solvent: float             # kg
    temperature: float              # K
    velocity: np.array=np.zeros(3)  # m/s
    position: np.array=np.zeros(3)  # m

    @staticmethod
    def from_mfs(solution, environment, radius, mass_fraction_solute, temperature,
                 velocity=np.zeros(3), position=np.zeros(3)):
        """Create a droplet from experimental conditions.

        Args:
            solution: parameters describing solvent+solute
            environment: parameters of gas surrounding droplet
            radius in metres
            mass_fraction_solute (MFS) (unitless)
            temperature in K
            velocity in metres/second (3-dimensional vector)
            position in metres (3-dimensional vector)
        """
        mass = 4*np.pi/3 * radius**3 * solution.density(mass_fraction_solute)
        mass_solvent = (1-mass_fraction_solute) * mass
        mass_solute = mass_fraction_solute * mass
        return UniformDroplet(solution, environment, mass_solute, mass_solvent, temperature, velocity, position)

    @property
    def state(self):
        """Variables that determine droplet state during its evolution, i.e. the independent variables."""
        return np.hstack((self.mass_solvent, self.temperature, self.velocity, self.position))

    @state.setter
    def state(self, new_state):
        """Update the state of the droplet by changing its independent variables.

        This is useful for reconstructing droplets at e.g. earlier points in time, or in its final
        equilibrium state.

        Args:
            new_state: the new set of independent variables.
        """
        try: self.mass_solvent, self.temperature, self.velocity, self.position = new_state
        except:
            x = new_state
            self.mass_solvent, self.temperature, self.velocity, self.position = x[0], x[1], x[2:5], x[5:]

    @property
    def complete_state(self):
        """All droplet variables, including both independent and dependent variables that completely
        determine all droplet properties.

        This form is ready for a row within a table specifying e.g. a droplet's trajectory.
        """
        return dict(mass = self.mass,
                    mass_solute = self.mass_solute,
                    mass_solvent = self.mass_solvent,
                    mass_fraction_solute = self.mass_fraction_solute,
                    mass_fraction_solvent = self.mass_fraction_solvent,
                    mole_fraction_solute = self.mole_fraction_solute,
                    mole_fraction_solvent = self.mole_fraction_solvent,
                    density = self.density,
                    radius = self.radius,
                    # This is too similar to radius so it can probably be omitted.
                    #diameter = self.diameter,
                    vapour_pressure = self.vapour_pressure,
                    temperature = self.temperature,
                    drag_coefficient = self.drag_coefficient,
                    reynolds_number = self.reynolds_number,
                    schmidt_number = self.schmidt_number,
                    prandtl_number = self.prandtl_number,
                    sherwood_number = self.sherwood_number,
                    nusselt_number = self.nusselt_number,
                    vx=self.velocity[0],
                    vy=self.velocity[1],
                    vz=self.velocity[2],
                    speed=self.speed,
                    x=self.position[0],
                    y=self.position[1],
                    z=self.position[2],
                    # I've assumed the user doesn't care about the derivatives, because these can be roughly inferred
                    # from the changes between frames. If more precise information is needed, then uncomment these.
                    #evaporation_rate = self.dmdt,
                    #dTdt = self.dTdt,
                    #ax = self.dvdt[0],
                    #ay = self.dvdt[1],
                    #az = self.dvdt[2]
                   )

    @property
    def mass(self):
        """Total mass of droplet (both solvent and solute components) in kg."""
        return self.mass_solute + self.mass_solvent

    @property
    def mass_fraction_solute(self):
        """Mass fraction of solute (i.e. the non-volatile component).
        NB: Should be zero for pure solvent."""
        return self.mass_solute / self.mass

    @property
    def mass_fraction_solvent(self):
        """Mass fraction of solvent."""
        return 1 - self.mass_fraction_solute

    @property
    def mole_fraction_solute(self):
        """Mole fraction of solute (i.e. the non-volatile component).
        NB: Should be zero for pure solvent."""
        return self.solution.mole_fraction_solute(self.mass_fraction_solute)

    @property
    def mole_fraction_solvent(self):
        """Mole fraction of solvent."""
        return 1 - self.mole_fraction_solute

    @property
    def density(self):
        """Droplet density in kg/m^3."""
        return self.solution.density(self.mass_fraction_solute)

    @property
    def radius(self):
        """Droplet radius in metres."""
        return (self.mass / (4*np.pi/3 * self.density))**(1/3)

    @property
    def diameter(self):
        """Droplet diameter in metres."""
        return 2*self.radius

    @property
    def vapour_pressure(self):
        """Vapour pressure at gas-liquid boundary in Pascals."""
        return self.solution.solvent_activity(self.mass_fraction_solute) * self.solution.solvent.equilibrium_vapour_pressure(self.temperature)

    @property
    def speed(self):
        """Magnitude of velocity vector in metres/second."""
        return np.linalg.norm(self.velocity)

    @property
    def relative_velocity(self):
        """Velocity relative to environment in metres/second."""
        return self.velocity - self.environment.velocity

    @property
    def relative_speed(self):
        """Magnitude of relative velocity vector in metres/second."""
        return np.linalg.norm(self.relative_velocity)

    @property
    def drag_coefficient(self):
        """Non-dimensional number describing strength of drag forces."""
        Re = self.reynolds_number
        if Re > 1000: return 0.424
        elif Re < 1e-12: return np.inf
        else: return (24 / Re) * (1 + Re**(2/3) / 6)

    @property
    def reynolds_number(self):
        """Non-dimensional number describing the type of fluid flow."""
        return self.environment.density * self.diameter * self.speed / self.environment.dynamic_viscosity

    @property
    def schmidt_number(self):
        """Non-dimensional number describing the ratio of momentum diffusivity to mass diffusivity."""
        D_function = self.solution.solvent.vapour_binary_diffusion_coefficient
        T_inf = self.environment.temperature
        D_inf = D_function(T_inf)
        return self.environment.dynamic_viscosity / (self.environment.density * D_inf)

    @property
    def prandtl_number(self):
        """Non-dimensional number describing the ratio of momentum diffusivity to thermal diffusivity."""
        return self.environment.specific_heat_capacity * self.environment.dynamic_viscosity / self.environment.thermal_conductivity

    @property
    def sherwood_number(self):
        """Non-dimensional number describing mass transfer."""
        Re = self.reynolds_number
        Sc = self.schmidt_number
        return 1 + 0.3 * Re**(1/2) * Sc**(1/3)

    @property
    def nusselt_number(self):
        """Non-dimensional number describing conductive heat transfer."""
        Re = self.reynolds_number
        Pr = self.prandtl_number
        return 1 + 0.3 * Re**(1/2) * Pr**(1/3)

    @property
    def dmdt(self):
        """Time derivative of mass, i.e. the rate of evaporation in kg/s."""
        Sh = self.sherwood_number

        D_function = self.solution.solvent.vapour_binary_diffusion_coefficient
        lam = D_function.lam

        T_inf = self.environment.temperature
        T = self.temperature
        D_inf = D_function(T_inf)

        # Apply temperature correction to diffusion coefficient appearing in mass flux.
        eps = 1e-8
        if np.abs(T_inf - T) < eps: C = 1 # ensure numerical stability as T -> T_inf
        else: C = (T_inf - T) / T_inf**(lam-1) * (2 - lam) / (T_inf**(2-lam) - T**(2-lam))
        D_eff = C * D_inf

        I = np.log((self.environment.pressure - self.vapour_pressure) /
                   (self.environment.pressure - self.environment.vapour_pressure))
        return 4*np.pi*self.radius * self.environment.density * D_eff * Sh * I

    @property
    def dTdt(self):
        """Time derivative of temperature from heat flux at the surface in K/s."""
        Nu = self.nusselt_number
        Gamma = stefan_boltzmann_constant

        r = self.radius
        m = self.mass
        rho = self.density

        T = self.temperature
        T_inf = self.environment.temperature
        K = self.environment.thermal_conductivity

        L = self.solution.solvent.specific_latent_heat_vaporisation(T)
        c = self.solution.solvent.specific_heat_capacity(T)
        r = self.radius

        return 3*K * (T_inf - T) * Nu / (c*rho*r**2) + L*self.dmdt / (c*m) - 3*Gamma * (T**4 - T_inf**4) / (c*rho*r)

    @property
    def dvdt(self):
        """Time derivative of velocity, i.e. its acceleration from Newton's second law in m/s^2."""
        rho_p = self.density
        rho_g = self.environment.density

        g = np.array((0, 0, -gravitational_acceleration))
        buoyancy = 1 - rho_g/rho_p
        acceleration = buoyancy*g

        C = self.drag_coefficient
        if np.isfinite(C):
            acceleration -= 3*C*rho_g * self.relative_speed * self.relative_velocity / (8*rho_p*self.radius)

        return acceleration

    @property
    def drdt(self):
        """Time derivative of droplet position, i.e. its velocity in m/s."""
        return self.velocity

    @property
    def dxdt(self):
        """Time derivative of all state variables as a vector."""
        return np.hstack((self.dmdt, self.dTdt, self.dvdt, self.drdt))

    @property
    def equilibrium_state(self):
        """Final state of the droplet once it has reached equilibrium."""

        temperature = self.environment.temperature
        velocity = np.zeros(3)
        position = np.zeros(3)

        # Find the right amount of solvent to equalise the vapour pressure across the gas-liquid boundary.
        vapour_pressure = self.environment.vapour_pressure
        solvent_activity = vapour_pressure / self.solution.solvent.equilibrium_vapour_pressure(temperature)
        mass_fraction_solute = self.solution.mass_fraction_solute_from_solvent_activity(solvent_activity)
        mass = self.mass_solute / mass_fraction_solute
        mass_solvent = mass - self.mass_solute

        return mass_solvent, temperature, velocity, position

    @property
    def equilibrium_droplet(self):
        """Final droplet once it has reached equilibrium."""
        return UniformDroplet(self.solution, self.environment, self.mass_solute, *self.equilibrium_state)

    def virtual_droplet(self, x):
        """Create droplet with new state variables.

        This is useful for creating droplets at future times.

        Args:
            x: new state variables (cf. Droplet.state function)
        """
        x = (x[0], x[1], x[2:5], x[5:])
        return UniformDroplet(self.solution, self.environment, self.mass_solute, *x)

    def copy(self):
        """Create an identical copy of this droplet."""
        return self.virtual_droplet(self.state.copy())

    def integrate(self, t, dt=1e-4, terminate_on_equilibration=False, eps=1e-4):
        """Integrate the droplet state forward in time.

        This solves an initial value problem with the current state as the initial conditions.
        The droplet state is updated to the final state after the integration.

        Args:
            t: total time to integrate over (s).
            dt: maximum timestep between frames in trajectory (s).
                NB: if numerical artifacts occur in the resulting trajectory, that suggests this parameter
                needs to be decreased.
            terminate_on_equilibration (default=False): if True, then the integration will stop if the
                evaporation rate falls below eps * the initial mass
            eps: threshold to use for the equilibration termination criterion.
        Returns:
            Trajectory of historical droplets showing how it reaches the new state.
        """
        from scipy.integrate import solve_ivp

        events = None
        if terminate_on_equilibration:
            m0 = self.mass
            equilibrated = lambda t,x: np.abs(self.virtual_droplet(x).dmdt) - eps*m0
            equilibrated.terminal = True
            events = [equilibrated]

        dxdt = lambda t,x: self.virtual_droplet(x).dxdt
        try:
            with np.errstate(divide='raise', invalid='raise'):
                trajectory = solve_ivp(dxdt, (0, t), self.state, max_step=dt, events=events)
        except:
            raise RuntimeError('an unknown error occurred during the simulation - try running again with a smaller timestep') from None

        self.state = trajectory.y[:,-1]
        return trajectory

    def complete_trajectory(self, trajectory):
        """Get the trajectory of all variables (including dependent ones) from a simulation (i.e.
        the output of UniformDroplet.integrate).

        Args:
            trajectory: the output of UniformDroplet.integrate, which gives the trajectory of independent
                        variables only.
        Returns:
            A pandas dataframe detailing the complete droplet history. 
        """

        variables = self.complete_state
        for label in variables:
            variables[label] = np.empty(trajectory.t.size)

        for i,state in enumerate(trajectory.y.T):
            earlier_droplet = self.virtual_droplet(state)
            earlier_state = earlier_droplet.complete_state
            for label,value in earlier_state.items():
                variables[label][i] = value

        variables['time'] = trajectory.t

        return pd.DataFrame(variables)


# # 3. Running the simulation

# ## 3.1. Define the Graphical User Interface (GUI) for specifying droplet simulations

# +
import ipywidgets as widgets
from IPython.display import display, HTML

from tkinter import Tk
from tkinter import filedialog

class DropletSimulationGUI:
    """This class manages a graphical user interface that allows the user to set
    droplet conditions and initiate simulations.

    Create and run this GUI by executing the following in an input cell:

    >>> gui = DropletSimulationGUI()
    >>> gui.display()
    """

    def __init__(self):
        """Constructor describes the various widgets used in the interface."""

        # The list of droplet trajectories from various simulations run by user.
        self.trajectories = []

        # Formatting options for widgets.
        dropdown_layout = widgets.Layout(width='100%')
        slide_layout = widgets.Layout(width='100%')
        label_layout = widgets.Layout(width='100%', display='flex', justify_content='center')
        centered_layout = widgets.Layout(width='100%', justify_content='center')
        coord_layout = widgets.Layout(width='100%')

        solution_grid_layout = widgets.Layout(grid_template_columns='repeat(2, 30%)')
                                             #justify_content='center')
        coord_grid_layout = widgets.Layout(grid_template_columns='repeat(3, 30%)',
                                           justify_items='center', align_items='center')
        time_grid_layout = widgets.Layout(grid_template_columns='20% 30% 50%',
                                          justify_items='center', align_items='center')
        tab_grid_layout = widgets.Layout(grid_template_columns='repeat(2, 50%)',
                                         #justify_items='center', # uncomment to horizontally align buttons in the center 
                                         align_items='center')

        # Section labels.
        self.droplet_label = widgets.Label('Droplet properties', layout=label_layout)
        self.initial_label = widgets.Label('Initial conditions', layout=label_layout)
        self.ambient_label = widgets.Label('Ambient conditions', layout=label_layout)
        self.simulation_label = widgets.Label('Simulation parameters', layout=label_layout)
        self.labels = [self.droplet_label, self.initial_label, self.ambient_label, self.simulation_label]
        # Hack to change the formatting of the section labels so they are more prominent.
        display(HTML("<style>.bolded { font-weight: bold; font-size: large }</style>"))
        for label in self.labels: label.add_class("bolded")

        ## Droplet properties.

        self.solution_dropdown = widgets.Dropdown(
            options=all_solutions.keys(),
            value=list(all_solutions.keys())[0],
            description='Solution:',
            disabled=False,
            layout=dropdown_layout
        )

        self.density_dropdown = widgets.Dropdown(
            options=['regular', 'volume additivity'],
            value='regular',
            description='Density fit:',
            disabled=False,
            layout=dropdown_layout
        )
        self.density_dropdown.observe(self.changed_density_fit, 'value')

        self.profile_dropdown = widgets.Dropdown(
            options=['uniform', 'radial'],
            value='uniform',
            description='Profile:',
            disabled=True,
            layout=dropdown_layout
        )
        self.profile_dropdown.observe(self.changed_profile, ['value', 'disabled'])

        self.npoints_select = widgets.BoundedIntText(
            value=100,
            min=10,
            max=10000,
            step=10,
            description='npoints:',
            disabled=True,
            layout=dropdown_layout
        )

        #self.solution_choices = widgets.HBox([self.solution_dropdown, self.density_dropdown])
        self.solution_choices = widgets.GridBox([self.solution_dropdown,
                                                 self.density_dropdown,
                                                 self.profile_dropdown,
                                                 self.npoints_select],
                                                 layout=solution_grid_layout)

        ## Initial conditions.

        self.initial_mfs_slider = widgets.FloatSlider(
            description='MFS', layout=slide_layout,
            value=0.05, min=0, max=1, step=0.01, readout_format='.2f')

        self.initial_R_slider = widgets.FloatSlider(
            description='R / µm', layout=slide_layout,
            value=15, min=0.1, max=500, step=0.1, readout_format='.1f')

        self.initial_T_slider = widgets.FloatSlider(
            description='T / ℃', layout=slide_layout,
            value=(body_temperature - T_freezing), min=0, max=100, step=0.1, readout_format='.1f')

        self.initial_x = widgets.FloatText(value=0, step=0.1, description='x / m', layout=coord_layout)
        self.initial_y = widgets.FloatText(value=0, step=0.1, description='y / m', layout=coord_layout)
        self.initial_z = widgets.FloatText(value=0, step=0.1, description='z / m', layout=coord_layout)
        self.initial_vx = widgets.FloatText(value=0, step=0.1, description='vx / m/s', layout=coord_layout)
        self.initial_vy = widgets.FloatText(value=0, step=0.1, description='vy / m/s', layout=coord_layout)
        self.initial_vz = widgets.FloatText(value=0, step=0.1, description='vz / m/s', layout=coord_layout)
        self.initial_coordinates = widgets.GridBox([self.initial_x, self.initial_y, self.initial_z,
                                                    self.initial_vx, self.initial_vy, self.initial_vz],
                                                    layout=coord_grid_layout)
        ## Ambient conditions.

        self.ambient_T_slider = widgets.FloatSlider(
            description='T / ℃', layout=slide_layout,
            value=(standard_temperature - T_freezing), min=0, max=100, step=0.1, readout_format='.1f')

        self.ambient_RH_slider = widgets.FloatSlider(
            description='RH', layout=slide_layout,
            value=0.5, min=0, max=1, step=0.01, readout_format='.2f')

        self.gas_vx = widgets.FloatText(value=0, step=0.1, description='vx / m/s', layout=coord_layout)
        self.gas_vy = widgets.FloatText(value=0, step=0.1, description='vy / m/s', layout=coord_layout)
        self.gas_vz = widgets.FloatText(value=0, step=0.1, description='vz / m/s', layout=coord_layout)
        self.gas_velocity = widgets.GridBox([self.gas_vx, self.gas_vy, self.gas_vz],
                                             layout=coord_grid_layout)

        ## Simulation parameters.

        self.time_selection = widgets.BoundedFloatText(
            value=10, min=0.1, max=1e6, step=0.1,
            description='time / s',
            layout=slide_layout
        )

        self.stop_checkbox = widgets.Checkbox(
            value=True,
            description='terminate on equilibration',
            indent=True,
            layout=slide_layout
        )
        self.stop_checkbox.observe(self.clicked_stop_checkbox, 'value')

        self.stop_threshold_slider = widgets.FloatLogSlider(
            description='threshold', layout=slide_layout,
            value=1e-3, base=10, min=-8, max=-1, step=1, readout_format='.1g'
        )

        #self.time_choices = widgets.HBox([self.time_selection, self.stop_checkbox, self.stop_threshold_slider])
        self.time_choices = widgets.GridBox([self.time_selection, self.stop_checkbox, self.stop_threshold_slider],
                                             layout=time_grid_layout)

        self.timestep_slider = widgets.FloatLogSlider(
            description='timestep / s', layout=slide_layout,
            value=1e-2, base=10, min=-5, max=-1, step=0.1, readout_format='.1g'
        )

        ## The start simulation button.

        self.run_button = widgets.Button(
            value=False,
            description='Simulate droplet',
            button_style='info',
            tooltip='Run simulation',
            icon='tint',
            style={'font_weight': 'bold'},
            layout=widgets.Layout(width='50ex')
        )

        self.run_button.on_click(self.run_simulation)

        self.centered_run_button = widgets.HBox([self.run_button], layout=centered_layout)

        ## Space to show warnings etc.

        self.status_space = widgets.HBox([], layout=centered_layout)

        ## Create tabs for showing plots of the various outputs.

        tab_labels = ['radius', 'mass', 'temperature', 'position']
        self.tabs = {label: widgets.GridBox([], layout=tab_grid_layout) for label in tab_labels}

        self.output_tabs = widgets.Tab(children=list(self.tabs.values()))
        for i, label in enumerate(tab_labels):
            self.output_tabs.set_title(i, label)

        ## Putting the entire interface together.

        self.display_widgets = [
                self.droplet_label,
                self.solution_choices,
                self.initial_label,
                self.initial_mfs_slider,
                self.initial_R_slider,
                self.initial_T_slider,
                self.initial_coordinates,
                self.ambient_label,
                self.ambient_T_slider,
                self.ambient_RH_slider,
                self.gas_velocity,
                self.simulation_label,
                self.time_choices,
                self.timestep_slider,
                self.centered_run_button,
                self.status_space,
                self.output_tabs
        ]

    def changed_density_fit(self, dropdown):
        """Callback when the 'density fit' dropdown is changed.

        If volume additivity is selected, then we can enable the options to allow the concentration
        profile to vary.

        Args:
            dropdown: the dropdown widget that is changed. This is not used, but it is needed
                      to assign this function as a callback to clicking the dropdown.
        """
        self.profile_dropdown.disabled = not self.density_dropdown.value == 'volume additivity'

        if not self.density_dropdown.value == 'volume additivity':
            self.profile_dropdown.value = 'uniform'

    def changed_profile(self, dropdown):
        """Callback when the 'profile' dropdown is changed.

        Unless we are doing inhomogeneous concentration profiles, we do not need to let the user
        specify the number of points so this function disables that widget.

        Args:
            dropdown: the dropdown widget that is changed. This is not used, but it is needed
                      to assign this function as a callback to clicking the dropdown.
        """
        self.npoints_select.disabled = self.profile_dropdown.value == 'uniform' or \
                                       not self.density_dropdown.value == 'volume additivity'

    def clicked_stop_checkbox(self, checkbox):
        """Callback when the 'terminate on equilibration' checkbox is clicked.

        We disable the termination threshold slider when this box is not clicked, because that parameter
        has no purpose outside of that termination condition.

        Args:
            checkbox: the checkbox widget that is clicked. This is not used, but it is needed
                      to assign this function as a callback to clicking the checkbox.
        """
        self.stop_threshold_slider.disabled = not self.stop_checkbox.value

    @property
    def gui_state(self):
        """Get the state of *all* GUI settings describing the parameters of a simulation.

        It is useful to keep these with trajectories for reference (so e.g. old simulations can be repeated).

        Return:
            dictionary containing all parameters defining a simulation.
        """
        return dict(
            # Droplet properties
            solution = self.solution_dropdown.value,
            density_fit = self.density_dropdown.value,
            profile = self.profile_dropdown.value,
            npoints = self.npoints_select.value,

            # Initial conditions of droplet.
            initial_radius = self.initial_R_slider.value,
            initial_mfs = self.initial_mfs_slider.value,
            initial_temperature = self.initial_T_slider.value,
            initial_velocity = np.array((self.initial_vx.value, self.initial_vy.value, self.initial_vz.value)),
            initial_position = np.array((self.initial_x.value, self.initial_y.value, self.initial_z.value)),

            # Environment (i.e. the gas phase).
            ambient_temperature = self.ambient_T_slider.value,
            ambient_RH = self.ambient_RH_slider.value,
            gas_velocity = np.array((self.gas_vx.value, self.gas_vy.value, self.gas_vz.value)),

            # Simulation parameters.
            time = self.time_selection.value,
            timestep = self.timestep_slider.value,
            terminate_on_equilibration = self.stop_checkbox.value,
            terminate_threshold = self.stop_threshold_slider.value
        )

    def run(self, solution, density_fit, profile, npoints,
            initial_radius, initial_mfs, initial_temperature, initial_velocity, initial_position,
            ambient_temperature, ambient_RH, gas_velocity,
            time, timestep, terminate_on_equilibration, terminate_threshold):
        """
        Run a new simulation.

        Args:
            solution: string specifying the solution i.e. (composition of droplet).
            density_fit: string specifying the type of density fit to use (valid values are 'regular' or
                         'volume additivity').
            profile: string either 'uniform' for uniform droplets or 'radial' for droplets with a
                     radially varying concentration profile. For the latter, the 'volume additivity'
                     density fit must be selected.
            npoints: number of points to use for radially varying concentration profiles (only used if
                     profile='radial').
            initial_radius: initial droplet radius in microns.
            initial_mfs: initial mass fraction of solute in droplet, should be between 0 and 1.
            initial_temperature: initial temperature of droplet in Celsius.
            initial_velocity: a 3-dimensional vector specifying the initial droplet velocity in metres/second.
            initial_position: a 3-dimensional vector specifying the initial droplet position in metres.
            ambient_temperature: temperature of surrounding gas in Celsius.
            ambient_RH: relative humidity of surrounding gas, should be between 0 and 1.
            gas_velocity: a 3-dimensional vector specifying the velocity of background gas in metres/second.
            time: the total time to run the simulation for in seconds.
            timestep: the timestep in the simulation in seconds. A smaller number will make the simulation
                      more accurate, but also take longer. NB: this is technically the *maximum* timestep,
                      as the integration algorithm will sometimes decide to take smaller timesteps.
            terminate_on_equilibration: a boolean; if True then the simulation will terminate early if
                      the droplet has stopped evolving (as determined by a threshold in the evaporation rate)
            terminate_threshold: the threshold for the previous termination condition. If the evaporation rate
                      divided by the initial mass falls below this threshold, then the simulation terminates.
        Returns:
            The final droplet state.
            The trajectory of (independent) droplet variables from which its complete history can be
                reconstructed.
        """

        ## Unit conversions.

        initial_radius *= 1e-6             # convert to metres
        initial_temperature += T_freezing  # convert to Kelvin
        ambient_temperature += T_freezing  # convert to Kelvin

        ## Set up droplet and environment.

        if density_fit == 'regular':
            solution = all_solutions[solution]
        elif density_fit == 'volume additivity':
            solution = all_solutions_volume_additivity[solution]
        else:
            raise ValueError('unrecognised density fit: %s' % density_fit)

        gas = Atmosphere(ambient_temperature, ambient_RH, velocity=gas_velocity)

        if profile == 'uniform':
            droplet = UniformDroplet.from_mfs(solution,
                                              gas,
                                              initial_radius,
                                              initial_mfs,
                                              initial_temperature,
                                              initial_velocity,
                                              initial_position)
        elif profile == 'radial':
            raise NotImplementedError('varying droplets not yet implemented!')
        else:
            raise ValueError('unrecognised concentration profile: %s' % profile)

        trajectory = droplet.integrate(time,
                                       timestep,
                                       terminate_on_equilibration,
                                       terminate_threshold)

        return droplet, trajectory

    def run_simulation(self, button):
        """Run a simulation with the currently selected settings.

        This function is called when the simulate button is clicked by the user.

        Args:
            button: the button widget that is clicked. This is not used, but it is needed
                    to assign this function as a callback to clicking the button.
        """

        ## Run the simulation with the current settings.
        settings = self.gui_state
        label_layout = widgets.Layout(width='100%', display='flex', justify_content='center')
        try:
            self.status_space.children = (widgets.Label("Running simulation...", layout=label_layout),)
            droplet, trajectory = self.run(**settings)
        except Exception as e:
            self.status_space.children = (widgets.Button(description=str(e),
                                                         disabled=True,
                                                         button_style='warning',
                                                         icon='exclamation-triangle',
                                                         style={'font_weight': 'bold'},
                                                         layout=label_layout),)
            return

        # On successful run, clear the status space.
        self.status_space.children = tuple()

        ## Save and output this new trajectory with some control widgets.

        # Add a new plot for this trajectory in each tab.
        plots = []
        for tab in self.tabs.values():
            plots += [widgets.Output()]
            tab.children = tab.children + (plots[-1],)
        self.plot_trajectory(droplet, trajectory, *plots)

        # Add a button so that the user can save this data to a file.
        save_button = widgets.Button(
            description='Save trajectory',
            button_style='success',
            tooltip='Save data to file',
            icon='download',
            style={'font_weight': 'bold'},
            layout=widgets.Layout(width='25ex')
        )
        save_button.on_click(self.save_trajectory)

        # Add a button so that the user can get more information about this trajectory.
        info_button = widgets.Button(
            description='Info',
            button_style='info',
            tooltip='More information about this trajectory (e.g. the initial conditions)',
            icon='info',
            style={'font_weight': 'bold'},
            layout=widgets.Layout(width='25ex')
        )
        info_button.on_click(self.info_trajectory)

        # Add a button so that the user can clear this data.
        remove_button = widgets.Button(
            description='Remove',
            button_style='danger',
            tooltip='Remove this trajectory',
            icon='times',
            style={'font_weight': 'bold'},
            layout=widgets.Layout(width='25ex')
        )
        remove_button.on_click(self.remove_trajectory)

        # Create a space for the buttons to the right of the plot.
        buttons = widgets.HBox([widgets.VBox([save_button, info_button, remove_button])],
                                layout=widgets.Layout(align_items='center'))
        for tab in self.tabs.values():
            tab.children = tab.children + (buttons,)

        # Save the data in case the user wants to dump it to a file (or otherwise use it).
        self.trajectories += [dict(droplet=droplet, trajectory=trajectory, settings=settings,
                                   plots=plots, buttons=buttons)]

    def save_trajectory(self, button):
        """Save a simulation to a file.

        This function is called when the save button is clicked by the user.

        Args:
            button: the save button widget that is clicked. We use its identity to determine which
                    trajectory to save.
        """

        # Find the trajectory in our saved list.
        for index,trajectory in enumerate(self.trajectories):
            if button in trajectory['buttons'].children[0].children:
                break
        else: # if the trajectory was not found, we have a problem!
            raise LookupError('could not find trajectory!')

        ## Prompt the user for a filepath to save the simulation, using Tkinter

        # We do not want to create a new window, so keep the Tkinter root window from appearing
        Tk().withdraw()
        # Prompt user for save path.
        filename = filedialog.asksaveasfilename(filetypes = (("comma separated values files", "*.csv"),
                                                             ("all files", "*.*")))

        ## Save the trajectory.

        if filename:
            droplet = trajectory['droplet']
            history = trajectory['trajectory']
            settings = trajectory['settings']
            settings = '\n'.join(['# %s: %s' % (key, value) for key,value in settings.items()])
            data = droplet.complete_trajectory(history)
            with open(filename, 'w') as f:
                f.write(settings)
                f.write('\n')
                data.to_csv(f, index=False)

    def info_trajectory(self, button):
        """Show more information about a previously run simulation.

        This function is called when the info button is clicked by the user.

        Args:
            button: the info button widget that is clicked. We use its identity to determine which
                    trajectory to show more information on.
        """

        # Find the trajectory in our saved list.
        for index,trajectory in enumerate(self.trajectories):
            if button in trajectory['buttons'].children[0].children:
                break
        else: # if the trajectory was not found, we have a problem!
            raise LookupError('could not find trajectory!')

        buttons = trajectory['buttons']
        settings = trajectory['settings']

        # If we're already showing the information, remove it.
        if len(buttons.children) > 1:
            buttons.children = tuple(list(buttons.children)[:-1])
        # Otherwise, display the information:
        else:
            css_style = '<style>p{word-wrap: break-word; line-height: 1.2em; margin-left: 1em}</style>'
            formatted_settings = '\n'.join(['<p>%s: %s</p>' % (key, value) for key,value in settings.items()])
            text = widgets.HTML(css_style + formatted_settings)
            buttons.children = buttons.children + (text,)

    def remove_trajectory(self, button):
        """Remove a simulation from the saved list.

        This function is called when the close button is clicked by the user.

        Args:
            button: the close button widget that is clicked. We use its identity to determine which
                    trajectory to remove.
        """

        # Find and remove the trajectory in our saved list.
        for index,trajectory in enumerate(self.trajectories):
            if button in trajectory['buttons'].children[0].children:
                break
        else: # if the trajectory was not found, we have a problem!
            raise LookupError('could not find trajectory!')
        self.trajectories.pop(index)

        # Remove the output rows/widgets from the display.
        for tab in self.tabs.values():
            # Remove any plots.
            tab.children = tuple(child for child in tab.children if child not in trajectory['plots'])
            # Remove the buttons.
            tab.children = tuple(child for child in tab.children if child is not trajectory['buttons'])

    def plot_trajectory(self, droplet, trajectory, out_radius, out_mass, out_temperature, out_position):
        """Plot summary data for a new trajectory.
        
        Args:
            droplet: the Droplet object, containing important information about the trajectory.
            trajectory: the trajectory of independent variables itself.
            out_radius: output pane for the radius trajectory plot.
            out_mass: output pane for the mass trajectory plot.
            out_temperature: output pane for the temperature trajectory plot.
            out_position: output pane for the position trajectory plot.
        """

        # Get the equilibrium state so we can compare how far/close to it we are.
        try:
            # If RH = 1, there is no equilibrium state, so we need to catch this.
            with np.errstate(divide='raise'):
                equilibrium = droplet.equilibrium_droplet
        except:
            # This should only happen if RH = 1.
            assert droplet.environment.relative_humidity == 1
            equilibrium = droplet.copy()
            equilibrium.mass_solvent = np.inf

        data = droplet.complete_trajectory(trajectory)

        t = data['time']
        radius = data['radius']
        mass_solvent = data['mass_solvent']
        temperature = data['temperature']
        x = data['x']
        z = data['z']

        figsize = (6, 4)

        with out_radius:
            plt.figure(figsize=figsize)
            plt.plot(t, radius, label='trajectory')
            plt.axhline(y=equilibrium.radius, ls='dashed', label='equilibrium')
            plt.legend(loc='best')
            plt.xlabel('t (s)')
            plt.ylabel('R (m)')
            plt.show()

        with out_mass:
            plt.figure(figsize=figsize)
            plt.plot(t, droplet.mass_solute + mass_solvent, label='total')
            plt.plot(t, mass_solvent, label='solvent')
            plt.axhline(y=equilibrium.mass, ls='dashed', label='equilibrium')
            plt.axhline(y=droplet.mass_solute, ls='dotted', label='solute')
            plt.legend(loc='best')
            plt.xlabel('t (s)')
            plt.ylabel('m (kg)')
            plt.show()

        with out_temperature:
            plt.figure(figsize=figsize)
            plt.plot(t, temperature - T_freezing, label='trajectory')
            plt.axhline(y=(equilibrium.temperature - T_freezing), ls='dashed', label='equilibrium')
            plt.legend(loc='best')
            plt.xlabel('t (s)')
            plt.ylabel('T (℃)')
            plt.show()

        with out_position:
            plt.figure(figsize=figsize)
            plt.plot(x, z)
            plt.xlabel('x (m)')
            plt.ylabel('z (m)')
            plt.show()

    def display(self):
        """Show the GUI display."""
        display(*self.display_widgets)

    def clear(self):
        """Clear all previous simulations."""
        self.trajectories = []
        for tab in self.tabs.values():
            tab.children = []


# -

# ## 3.2. Executing the graphical program to run simulations

gui = DropletSimulationGUI()
gui.display()

# # Appendix: Differential Equations

# $\renewcommand\vec[1]{\mathbf{#1}}$
# The equations of motion in Xie et al (2007) read:
#
# \begin{align}
#   \frac{d r_p}{dt} =& \;
#   \mathrm{Sh}
#   \frac{C D_\infty}{\rho_p r_p} \frac{M_v p}{R T_\infty}
#   \ln{\left( \frac{p - p_{va}}{p - p_{v^\infty}} \right)}
#   = f_1(r_p, T_p, \vec{V}_p),
#   \\
#   \frac{d T_p}{dt} =& \;
#   4\pi r_p^2
#   \, K_g \frac{T_\infty - T_p}{c_p m_p r_p} \mathrm{Nu}
#   - \frac{L_v I}{c_p m_p}
#   - 4\pi r_p^2 \frac{\Gamma (T_p^4 - T_\infty^4)}{c_p m_p}
#   \\ =& \;
#   3 K_g \frac{T_\infty - T_p}{c_p \rho_p r_p^2} \mathrm{Nu}
#   + \frac{L_v \dot{m}_p}{c_p m_p}
#   - \frac{3 \Gamma (T_p^4 - T_\infty^4)}{c_p \rho_p r_p}
#   = f_2(r_p, T_p, \vec{V}_p),
#   \\
#   \frac{d \vec{V}_p}{dt} =& \;
#   \vec{g} \left( 1 - \frac{\rho_p}{\rho_g} \right)
#   - \frac{3 C_d \rho_g |\vec{V}_p - \vec{V}_g| (\vec{V}_p - \vec{V}_g)}{8 \rho_p r_p}
#   = f_3(r_p, T_p, \vec{V}_p),
#   \\
#   \frac{d \vec{x}_p}{dt} =& \;
#   \vec{V}_p
#   = f_4(\vec{V}_p).
# \end{align}
#
# These equations implicitly assume a further relationship to solve these.
# To see this, we divide both the numerator/denominator inside the logarithm in the first equation through by $p$ to obtain
#
# \begin{equation}
#   \frac{p - p_{va}}{p - p_{v^\infty}}
#   =
#   \frac{p (1 - \chi_w(r_p^+))}{p - p_{v^\infty}}
#   =
#   \frac{1 - \chi_w(r_p^+)}{1 - \chi_w(\infty)},
# \end{equation}
#
# which requires a further equation for $\chi_w(r_p^+)$ (the water/solvent mole fraction) for closure.
# Xie et al (2007) assume Raoult's law, i.e.
#
# \begin{equation}
#   p_{va} = \chi_w(r_p^-) \, p_{va,s}(T_p).
# \end{equation}
#
# Mechanical balance ensures that $p_{va}$ is continous across the liquid-vapour phase boundary, and so $p_{va} = \chi_w(r_p^+) \, p$ using the ideal gas law.
#
# The equilibrium vapour pressure for the solvent at the phase boundary can be calculated using the Clausius-Clapeyron relation, leading to Buck's equation for $p_{va}(T_p)$.

# #### Constant density assumption?
#
# The first equation implicitly assumes that the density is constant.
# To see this, we rewrite this equation as the following two subequations:
#
# \begin{align}
#   \frac{d r_p}{dt} =& \;
#   \frac{1}{4\pi r_p^2 \rho_p} \frac{d m_p}{dt},
#   \\
#   \frac{d m_p}{dt} =& \;
#   \mathrm{Sh} \,
#   4\pi r_p \, C D_\infty \, \frac{M_v p}{R T_\infty}
#   \ln{\left( \frac{p - p_{va}}{p - p_{v^\infty}} \right)}.
# \end{align}
#
# Next, consider the expression for total droplet mass:
#
# \begin{equation}
#   m_p = 4\pi \int_0^{r_p} \rho_p(r') r'^2 \, dr'.
# \end{equation}
#
# When the density is uniform inside the droplet $\rho_p(r') = \rho_p$ then we obtain
#
# \begin{equation}
#   m_p = \frac{4\pi r_p^3}{3} \rho_p,
# \end{equation}
#
# and so the time-derivative yields
#
# \begin{equation}
#   \frac{d m_p}{dt}
#   =
#   4\pi r_p^2 \rho_p \frac{d r_p}{dt}
#   +
#   \frac{4\pi r_p^3}{3} \frac{d \rho_p}{dt}.
# \end{equation}
#
# If we make the density constant in time then this yields the expression for $\dot{r}_p$ above.
# It is possible this is approximately valid if the density only varies slowly.
# More generally, if we have a varying density profile then we obtain
#
# \begin{equation}
#   \frac{d m_p}{dt}
#   =
#   4\pi r_p^2 \rho_p(r_p) \frac{d r_p}{dt}
#   +
#   4\pi \int_0^{r_p} \frac{\partial \rho_p(r')}{\partial t} r'^2 \, dr'.
# \end{equation}
