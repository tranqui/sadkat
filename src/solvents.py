# # 2. Physical parameters entering model

# ## 2.1 Solvents

# ### 2.1.1. Data structures for general solvents.

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

# ### 2.1.2. Properties of pure water

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

