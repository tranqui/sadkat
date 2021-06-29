# + ignore="True"
from preamble import *
# -

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
    surface_tension : object                     # N/m

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

# ### 2.1.2 Kelvin Effect


# +
def surface_tension(A, B, C, D, E, T_crit, T):
    """Takes fitting parameters and tempertature (K) and returns surface tension (N/m).
    http://ddbonline.ddbst.de/DIPPR106SFTCalculation/DIPPR106SFTCalculationCGI.exe
    """

    T_r = T / T_crit

    power = B + (C * T_r) + (D * T_r ** 2) + (E * T_r ** 3)

    sigma = (A * (1 - T_r) ** power) / 1000 # convert from mN/m to N/m


    return sigma

# cf. https://www.e-education.psu.edu/meteo300/node/676
def kelvin_effect(solvent_surface_tension, solvent_density, solvent_molar_mass, T, P_vap_flat, droplet_radius):
    """Takes solvent durface tension, density, molar mass, temperature, equlibrium vapour pressure and droplet radius
    and returns vapour pressure of curved surface.
    """

    n_L = solvent_density / solvent_molar_mass

    P_vap_curved =  P_vap_flat * np.exp( (2 * solvent_surface_tension) / (n_L * gas_constant * T * droplet_radius))

    return P_vap_curved
# -


# Plot kelvin effect with approx value of surface tension and p_vap_flat and compare to reference: https://www.e-education.psu.edu/meteo300/node/676
# This shows that it is only really significant for very small droplets, but I think its valuable and relatively easy to include.

# +
R_range = np.arange(1e-9,30e-9, 1e-10)

curve_range = kelvin_effect(0.073, 997, 0.018, 293, 2300, R_range )

plt.plot(R_range/1e-9,curve_range/2300)
plt.xlabel('radius / nm')
plt.ylabel('P_vap / P_vap_flat')
# -

# ### 2.1.3. Properties of pure water

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

# ? not sure where this is from
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

    return P


def surface_tension_water(T):
    """ parameters from     http://ddbonline.ddbst.de/DIPPR106SFTCalculation/DIPPR106SFTCalculationCGI.exe
        Tc  	 Tmin 	 Tmax
        647.3 	 233 	 643
        A,B,C,D,E = 134.15,1.6146,-2.035,1.5598,0
    """

    return surface_tension(134.15,1.6146,-2.035,1.5598,0, 647.3, T)


Water = Solvent(molar_mass_water,
                density_water,
                specific_heat_capacity_water,
                specific_latent_heat_water,
                equilibrium_vapour_pressure_water,
                VapourBinaryDiffusionCoefficient(0.2190e-4, T_freezing, 1.81),
                surface_tension_water)
# -

# Sanity check the water properties by plotting them below:

# +
if __name__ == '__main__':
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
    ax3.set_ylabel('equilibrium\nvapour pressure (Pa)')

    ax4.plot(T_C, Water.vapour_binary_diffusion_coefficient(T_C+T_freezing))
    ax4.set_xlabel('T (℃)')
    ax4.set_ylabel('$D_\infty$ (m$^2$/s)')

    plt.show()
# -

# +
plt.plot(T_C ,Water.surface_tension(T_C + T_freezing))
plt.xlabel('T (℃)')
plt.ylabel('Surface tension / N / m')

# -

# +
def rackett_equation(A, B, C, D, T):
    """ Rackett equation, takes coefficients and T (in Kelvin) and returns density.
    Taken from http://ddbonline.ddbst.de/DIPPR105DensityCalculation/DIPPR105CalculationCGI.exe?component=Methanol
    """

    density = (A) / (B ** (1 + ((1 - (T / C))) ** D))

    return density
# -

# +
def antoine_equation(T, A, B, C):
    """T in deg C, returns P in whatever parameters are for"""

    P =  10 ** ( A - (B / (T + C)))

    return P
# -

# +
# NIST latent heat of vapourisation parameterisation
#
# ![image.png](attachment:image.png)
# -

#+
def enthalpy_vapourisation(A, alpha, beta, T_c, T):
    """Rreceived T in K and coefficients. Returns H_vap in J/mol"""
    T_r = T / T_c

    H_vap = 1000 * ( A * np.exp( - alpha * T_r) * (1 - T_r) ** beta )

    return H_vap
# -

# +
def temp_K_to_C(temperature_K):
    return temperature_K - 273.15
# -

# ### 2.1.4. Properties of pure ethannol

# +
molar_mass_ethanol = 46.07 # g/mol

def density_ethanol(temperature):
    """ densisty as function of temperature
    http://ddbonline.ddbst.de/DIPPR105DensityCalculation/DIPPR105CalculationCGI.exe?component=2-ethanol
    T_min, t_max = 191,513

    Args:
        temperature: temperature in Kelvin.
    Returns:
        The ethanol density in kg/m^3.
    """

    A, B, C, D = 99.3974, 0.310729, 513.18, 0.30514


    density = rackett_equation(A, B, C, D, temperature)

    return density

def specific_heat_capacity_ethanol(T):
    # https://webbook.nist.gov/cgi/cbook.cgi?ID=C64175&Mask=2
    #using constant value and converting from molar to specific
    molar_mass_ethanol = 46.07 # g/mol
    molar_heat_capacity_ethanol = 112.3 # J/mol/K

    specific_heat_capacity = molar_heat_capacity_ethanol / molar_mass_ethanol * 1000

    return specific_heat_capacity

def specific_latent_heat_ethanol(T):
    """Temperature (K)	298. - 469 Temperature limits.
        A (kJ/mol)	50.43
        α	-0.4475
        β	0.4989
        Tc (K)	513.9
        Reference	Majer and Svoboda, 1985
        https://webbook.nist.gov/cgi/cbook.cgi?ID=C64175&Mask=4
    """

    A, alpha, beta, T_c = 50.43, -0.4475, 0.4989, 513.9
    molar_mass_ethanol = 46.07 # g/mol

    H_vap = enthalpy_vapourisation(50.43, -0.4475, 0.4989, 513.9, T) # J / Kg

    H_vap = H_vap / molar_mass_ethanol

    return H_vap

def equilibrium_vapour_pressure_ethanol(temperature_K):
    """Take temperature in kelvin and returns vapour pressure in Pa, probabaly a bit inaccurate below 20 deg c"""
    # https://webbook.nist.gov/cgi/inchi?ID=C71238&Mask=4&Type=ANTOINE&Plot=on

    P_vap = antoine_equation(temperature_K,5.37229, 1670.409, -40.191) * 100000 # coeffs for bar so have converted


    return P_vap

def surface_tension_ethanol(T):
    """http://ddbonline.ddbst.de/DIPPR106SFTCalculation/DIPPR106SFTCalculationCGI.exe
                A 	         B 	         C 	           D 	     E 	     Tc 	 Tmin 	 Tmax
                131.38 	 5.5437 	 -8.4826 	 4.3164 	 0 	 516.2 	 180 	 513
    """
    return surface_tension(131.38, 5.5437, -8.4826, 4.3164,0,516.2,T)

Ethanol = Solvent(molar_mass_ethanol,
                density_ethanol,
                specific_heat_capacity_ethanol,
                specific_latent_heat_ethanol,
                equilibrium_vapour_pressure_ethanol,
                VapourBinaryDiffusionCoefficient(0.2190e-4, T_freezing, 1.81),
                surface_tension_ethanol)
# -

# +
T_C = np.linspace(0, 78, 79)
fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(ncols=2, nrows=2)

colour = 'k'

ax1.plot(T_C, Ethanol.density(T_C + T_freezing), c = colour)
ax1.set_xlabel('T (℃)')
ax1.set_ylabel('density (kg/m3)')

ax2.plot(T_C, Ethanol.specific_latent_heat_vaporisation(T_C + T_freezing), c = colour)
#ax2.axhline(y=2264.705, ls='dashed')
ax2.set_xlabel('T (℃)')
ax2.set_ylabel('L (kJ/kg)')

ax3.plot(T_C, Ethanol.equilibrium_vapour_pressure(T_C+T_freezing), c = colour)
ax3.set_xlabel('T (℃)')
ax3.set_ylabel('equilibrium vapour pressure (Pa)')

ax4.plot(T_C, Ethanol.vapour_binary_diffusion_coefficient(T_C+T_freezing), c = colour)
ax4.set_xlabel('T (℃)')
ax4.set_ylabel('$D_\infty$ (m$^2$/s)')

plt.show()
# -

# +
plt.plot(T_C ,Ethanol.surface_tension(T_C + T_freezing), c = 'k')
plt.xlabel('T (℃)')
plt.ylabel('Surface tension / N / m')
# -

# ### 2.1.5. Properties of pure propanol

# +
molar_mass_propanol = 60.0952 # g/mol

def density_propanol(temperature):
    """ densisty as function of temperature
    http://ddbonline.ddbst.de/DIPPR105DensityCalculation/DIPPR105CalculationCGI.exe?component=2-propanol
    T_min, t_max = 186,  507

    Args:
        temperature: temperature in Kelvin.
    Returns:
        The propanol density in kg/m^3.
    """


    A, B, C, D = 74.5237, 0.27342, 508.3, 0.235299


    density = rackett_equation(A, B, C, D, temperature)

    return density


def specific_heat_capacity_propanol(T):
    # https://webbook.nist.gov/cgi/cbook.cgi?ID=C71238&Mask=2
    #using constant value and converting from molar to specific

    molar_mass_propanol = 60.0952 # g/mol
    molar_heat_capacity_propanol = 145 # J/mol/K

    specific_heat_capacity = molar_heat_capacity_propanol / molar_mass_propanol * 1000

    return specific_heat_capacity

def specific_latent_heat_propanol(T):
    """Temperature (K)	298. - 390.
        A (kJ/mol)	52.06
        α	-0.8386
        β	0.6888
        Tc (K)	536.7
        Reference	Majer and Svoboda, 1985

        https://webbook.nist.gov/cgi/cbook.cgi?ID=C64175&Mask=4
        """

    A, alpha, beta, T_c = 52.06, -0.8386, 0.6888, 536.7
    molar_mass_propanol = 60.0952 # g/mol

    H_vap = enthalpy_vapourisation(50.43, -0.4475, 0.4989, 513.9, T) # J / Kg

    H_vap = H_vap / molar_mass_propanol

    return H_vap

def equilibrium_vapour_pressure_propanol(temperature_K):
    """Take temperature in kelvin and returns vapour pressure in Pa, probabaly a bit inaccurate below 20 deg c """
    # https://webbook.nist.gov/cgi/inchi?ID=C71238&Mask=4&Type=ANTOINE&Plot=on

    temperature_C = temperature_K - 273.15


    P_vap = antoine_equation(temperature_K,5.31, 1690.86, -51.80) * 100000 # coeffs for bar so have converted

    return P_vap

def surface_tension_propanol(T):
    """ http://ddbonline.ddbst.de/DIPPR106SFTCalculation/DIPPR106SFTCalculationCGI.exe
    DIPPR106 Equation Parameters (Surface Tension in mN/m, T in K)
    No.	 A 	 B 	 C 	 D 	 E 	 Tc 	 Tmin 	 Tmax
    (1) 	 46.507 	 0.90053 	 0 	 0 	 0 	 508.3 	 287 	 353
    """

    return surface_tension(46.507,0.90053,0,0,0,508.3, T)

Propanol = Solvent(molar_mass_propanol,
                density_propanol,
                specific_heat_capacity_propanol,
                specific_latent_heat_propanol,
                equilibrium_vapour_pressure_propanol,
                VapourBinaryDiffusionCoefficient(0.2190e-4, T_freezing, 1.81),
                surface_tension_propanol)
# -

# +
T_C = np.linspace(0, 100, 101)
fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(ncols=2, nrows=2)

colour = 'r'

ax1.plot(T_C, Propanol.density(T_C + T_freezing), c = colour)
ax1.set_xlabel('T (℃)')
ax1.set_ylabel('density (kg/m3)')

ax2.plot(T_C, Propanol.specific_latent_heat_vaporisation(T_C + T_freezing), c = colour)
#ax2.axhline(y=2264.705, ls='dashed')
ax2.set_xlabel('T (℃)')
ax2.set_ylabel('L (kJ/kg)')

ax3.plot(T_C, Propanol.equilibrium_vapour_pressure(T_C+T_freezing), c = colour)
ax3.set_xlabel('T (℃)')
ax3.set_ylabel('equilibrium vapour pressure (Pa)')

ax4.plot(T_C, Propanol.vapour_binary_diffusion_coefficient(T_C+T_freezing), c = colour)
ax4.set_xlabel('T (℃)')
ax4.set_ylabel('$D_\infty$ (m$^2$/s)')

plt.show()
# -

# +
plt.plot(T_C ,Propanol.surface_tension(T_C + T_freezing), c = 'r')
plt.xlabel('T (℃)')
plt.ylabel('Surface tension / N / m')
# -

# ### 2.1.5. Properties of pure Butanol

# +
molar_mass_butanol = 74.12 # g/mol

def density_butanol(temperature):
    """ densisty as function of temperature
    http://ddbonline.ddbst.de/DIPPR105DensityCalculation/DIPPR105CalculationCGI.exe
    T_min, t_max = 213, 558

    Args:
        temperature: temperature in Kelvin.
    Returns:
        The butanol density in kg/m^3.
    """

    A, B, C, D =  9.87035, 0.0998069, 568.017, 0.126276


    density = rackett_equation(A, B, C, D, temperature)

    return density

def specific_heat_capacity_butanol(T):
    # https://webbook.nist.gov/cgi/cbook.cgi?ID=C71238&Mask=2
    #using constant value and converting from molar to specific

    molar_mass_butanol = 74.12 # g/mol
    molar_heat_capacity_butanol = 177 # J/mol/K

    specific_heat_capacity = molar_heat_capacity_butanol / molar_mass_butanol * 1000

    return specific_heat_capacity

def specific_latent_heat_butanol(T):
    """Temperature (K)	298. - 372.
        A (kJ/mol)	52.6
        α	-1.462
        β	1.0701
        Tc (K)	536.
        Reference	Majer and Svoboda, 1985
        https://webbook.nist.gov/cgi/cbook.cgi?ID=C78922&Mask=4
    """

    A, alpha, beta, T_c = 52.6, -1.462, 1.0701, 536
    molar_mass_butanol = 74.12 # g/mol

    H_vap = enthalpy_vapourisation(50.43, -0.4475, 0.4989, 513.9, T) # J / Kg

    H_vap = H_vap / molar_mass_butanol

    return H_vap

def equilibrium_vapour_pressure_butanol(temperature_K):
    """Take temperature in kelvin and returns vapour pressure in Pa, probabaly a bit inaccurate below 20 deg c """
    # https://webbook.nist.gov/cgi/cbook.cgi?ID=C71363&Mask=4&Type=ANTOINE&Plot=on

    temperature_C = temperature_K - 273.15

    P_vap = antoine_equation(temperature_K,4.55, 1351.56, -93.34) * 100000 # coeffs for bar so have converted

    return P_vap

def surface_tension_butanol(T):
    """
    DIPPR106 Equation Parameters (Surface Tension in mN/m, T in K)
    No.	 A 	 B 	 C 	 D 	 E 	 Tc 	 Tmin 	 Tmax
    (1) 	 72.697 	 3.0297 	 -4.2681 	 2.4776 	 0 	 562.9 	 238 	 543
    """
    return surface_tension(72.697, 3.0297, -4.2681, 2.4776,0, 562.9, T)

Butanol = Solvent(molar_mass_butanol,
                density_butanol,
                specific_heat_capacity_butanol,
                specific_latent_heat_butanol,
                equilibrium_vapour_pressure_butanol,
                VapourBinaryDiffusionCoefficient(0.2190e-4, T_freezing, 1.81),
                surface_tension_butanol)
# -

# +
T_C = np.linspace(0, 117, 118)
fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(ncols=2, nrows=2)

colour = 'g'

ax1.plot(T_C, Butanol.density(T_C + T_freezing), c = colour)
ax1.set_xlabel('T (℃)')
ax1.set_ylabel('density (kg/m3)')

ax2.plot(T_C, Butanol.specific_latent_heat_vaporisation(T_C + T_freezing), c = colour)
#ax2.axhline(y=2264.705, ls='dashed')
ax2.set_xlabel('T (℃)')
ax2.set_ylabel('L (kJ/kg)')

ax3.plot(T_C, Butanol.equilibrium_vapour_pressure(T_C+T_freezing), c = colour)
ax3.set_xlabel('T (℃)')
ax3.set_ylabel('equilibrium vapour pressure (Pa)')

ax4.plot(T_C, Butanol.vapour_binary_diffusion_coefficient(T_C+T_freezing), c = colour)
ax4.set_xlabel('T (℃)')
ax4.set_ylabel('$D_\infty$ (m$^2$/s)')

plt.show()
# -

# +
plt.plot(T_C ,Butanol.surface_tension(T_C + T_freezing), c = 'g')
plt.xlabel('T (℃)')
plt.ylabel('Surface tension / N / m')
# -
