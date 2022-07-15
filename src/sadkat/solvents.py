# + ignore="True"
from sadkat.common import *
# -

# # 2. Physical parameters entering model

# ## 2.1. Solvents

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
    """Standard fitting function to describe the temperature dependence
    of the binary diffusion coeffient for evaporated solvent in air.

    The functional form of this is taken from Xie et al., Indoor Air (2007).
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

# ### 2.1.2. Kelvin effect

# +
def surface_tension(A, B, C, D, E, T_crit, T):
    """Surface tension fitting function (as function of temperature).

    Source:
    http://ddbonline.ddbst.de/DIPPR106SFTCalculation/DIPPR106SFTCalculationCGI.exe

    Args:
        A, B, C, D, E: parameters of fit.
        T_crit: temperature of critical point in Kelvin.
        T: temperature(s) to evaluate surface tension at (Kelvin).
    Returns:
        Surface tension in N/m.
    """

    T_r = T / T_crit
    power = B + (C * T_r) + (D * T_r ** 2) + (E * T_r ** 3)
    sigma = (A * (1 - T_r) ** power) / 1000 # convert from mN/m to N/m

    return sigma

def kelvin_effect(solvent_surface_tension, solvent_density, solvent_molar_mass, T, P_vap_flat, droplet_radius):
    """Drop in vapour pressure due to surface curvature.

    Read more: https://www.e-education.psu.edu/meteo300/node/676

    Args:
        solvent_surface_tension: surface tension of droplet/medium interface (N/m).
        solvent_density: density of droplet (kg/m^3).
        solvent_molar_mass: molar mass of droplet (g/mol).
        T: temperature of droplet (K)
        P_vap_flat: vapour pressure of a flat interface under equivalent conditions (Pa).
        droplet_radius: radius of droplet (m).
    Returns:
        Vapour pressure of the curved interface (Pa).
    """

    n_L = solvent_density / solvent_molar_mass
    return P_vap_flat * np.exp( (2 * solvent_surface_tension) / (n_L * gas_constant * T * droplet_radius))
# -

# Plot kelvin effect with approx value of surface tension and p_vap_flat and compare to reference: https://www.e-education.psu.edu/meteo300/node/676
# This shows that it is only really significant for very small droplets, but I think its valuable and relatively easy to include.

if __name__ == '__main__':
    R_range = np.arange(1e-9, 30e-9, 1e-10)

    curve_range = kelvin_effect(0.073, 997, 0.018, 293, 2300, R_range)

    plt.figure()
    plt.plot(R_range/1e-9, curve_range/2300)
    plt.xlabel('radius / nm')
    plt.ylabel('P_vap / P_vap_flat')

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

#IAPWS-95 https://chemicals.readthedocs.io/chemicals.iapws.html
specific_heat_capacity_water = np.vectorize(lambda T: chemicals.iapws.iapws95_properties(T, standard_atmospheric_pressure)[5]) # J/kg/K

# Su, PCCP (2018)
specific_latent_heat_water = lambda T: 3.14566e6 - 2361.64 * T # J/kg

def equilibrium_vapour_pressure_water(T):
    """Using the Buck equation (from integrating the Clausius-Clapeyron equation).

    Args:
        T: temperature in Kelvin.
    Returns:
        The vapour pressure in Pascals.
    """
    T_C = T - T_freezing # Celsius
    return 1e3*0.61161 * np.exp((18.678 - (T_C / 234.5)) * (T_C / (257.14 + T_C))) # Pa

def surface_tension_water(T):
    """Surface tension of water as a function of temperature.

    Parameterisation of water surface tension from:
    http://ddbonline.ddbst.de/DIPPR106SFTCalculation/DIPPR106SFTCalculationCGI.exe

    Parameters of fit from source:
        Tc           Tmin    Tmax
        647.3        233     643
        A,B,C,D,E = 134.15, 1.6146, -2.035, 1.5598, 0

    Args:
        T: temperature in Kelvin.
    Returns:
        Surface temperature in N/m.
    """
    return surface_tension(134.15, 1.6146, -2.035, 1.5598, 0, 647.3, T)

Water = Solvent(molar_mass_water,
                density_water,
                specific_heat_capacity_water,
                specific_latent_heat_water,
                equilibrium_vapour_pressure_water,
                VapourBinaryDiffusionCoefficient(0.2190e-4, T_freezing, 1.81),
                surface_tension_water)
# -

# Sanity check the water properties by plotting them below:

if __name__ == '__main__':
    def plot_solvent_properties(solvent, T = np.linspace(T_freezing, T_freezing + 99.9, 201)):
        fig, axes = plt.subplots(ncols=2, nrows=3, sharex=True)
        plt.subplots_adjust(hspace=0)

        axes[0][0].plot(T, solvent.density(T), lw = 3, color = 'k')
        axes[0][0].set_ylabel(r'$\rho$ / kgm$^{-3}$')

        axes[0][1].plot(T, solvent.specific_heat_capacity(T) / 1e3, lw = 3, color = 'k')
        axes[0][1].set_ylabel('C$_p$ / kJK$^{-1}$kg$^{-1}$')

        axes[1][0].plot(T, solvent.specific_latent_heat_vaporisation(T) / 1e3, lw = 3, color = 'k')
        axes[1][0].set_ylabel('$\Delta$H$_{vaporisation}$ / kJkg$^{-1}$')

        axes[1][1].plot(T, solvent.equilibrium_vapour_pressure(T), lw = 3, color = 'k')
        axes[1][1].set_ylabel('P$_0$ / Pa')

        axes[2][0].plot(T, solvent.vapour_binary_diffusion_coefficient(T) / 1e-6, lw = 3, color = 'k')
        axes[2][0].set_ylabel(r'D$_{water, air}$ / 10$^{-6}$ $\times$ m$^2$s$^{-1}$')

        axes[2][1].plot(T, solvent.surface_tension(T) / 1e-3, lw = 3, color = 'k')
        axes[2][1].set_ylabel(r'$\sigma$ / mNm$^{-1}$')

    plot_solvent_properties(Water)

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

    plt.figure()
    plt.plot(T_C, Water.surface_tension(T_C + T_freezing))
    plt.xlabel('T (℃)')
    plt.ylabel('Surface tension / N / m');

# ### 2.1.4. Alcohols
# #### 2.1.4.1. Fit functions

# +
def rackett_equation(A, B, C, D, T):
    """Rackett equation used to parameterise density (vs temperature) for various solvents.

    Source:
    http://ddbonline.ddbst.de/DIPPR105DensityCalculation/DIPPR105CalculationCGI.exe

    Args:
        A, B, C, D: fit parameters.
        T: temperature in Kelvin.
    Returns:
        Density in kg/m^3
    """
    return A / (B ** (1 + ((1 - (T / C))) ** D))

def antoine_equation(T, A, B, C):
    """Semi-empirical fit function describing relation between vapour pressure and temperature.

    Args:
        T: temperature in Kelvin.
        A, B, C: fit parameters.
    Returns:
        Vapour pressure in Pa.
    """
    return 10 ** ( A - (B / (T + C)))

def enthalpy_vapourisation(A, alpha, beta, T_c, T):
    """Fit function used by NIST for enthalpy of vaporisation H_vap (at saturation pressure).

    Fit function is that described in Majer and Svoboda (1985), referenced to by NIST here:
    https://webbook.nist.gov/cgi/cbook.cgi?ID=C64175&Mask=4

    Args:
        A, alpha, beta: fit parameters of function.
        T_c: critical temperature (Kelvin).
        T: temperature(s) to evaluate H_vap at (Kelvin).
    Returns:
        Enthalpy of vapourisation H_vap in J/mol.
    """
    T_r = T / T_c
    return 1000 * ( A * np.exp( - alpha * T_r) * (1 - T_r) ** beta )
# -

# #### 2.1.4.2. Alcohol diffusion coefficients
#
# Reference diffusion values: https://www.sciencedirect.com/science/article/pii/B978032328658900010X

# +
def yaws_diffusion_polynomial(T, A, B, C):
    return 1e-4 * (A + B*T + C*T**2)

def vapour_binary_diffusion_coefficient_func(T, D_ref, T_ref, lam):
    return D_ref * (T / T_ref)**lam

yaws_alcohol_properties = [
    [40, 'CH4O', 'methyl_alcohol', 67_56_1, -0.12691, 7.2728E-04, 7.0516E-07, 200, 1500, 1,2, 0.0469, 0.1526, 2.5506],
    [175, 'C2H6O', 'ethyl_alcohol', 64_17_5, -0.10107, 5.6275E-04, 5.8314E-07, 200, 1500, 1,2, 0.0349, 0.1186, 2.0551],
    [353, 'C3H8O', 'propyl_alcohol', 71_23_8, -0.08697, 4.7538E-04, 5.0525E-07, 200, 1500, 1,2, 0.0284, 0.0997, 1.7629],
    [354, 'C3H8O', 'isopropyl_alcohol', 67_63_0, -0.08490, 4.7550E-04, 5.0399E-07, 200, 1500, 2, 0.0305, 0.1017, 1.7623],
    [689, 'C4H10O', 'butanol', 71_36_3, -0.07689, 4.1316E-04, 4.5143E-07, 200, 1500, 1,2, 0.0239, 0.0864, 1.5586],
    [690, 'C4H10O', 'isobutanol', 78_83_1, -0.07577, 4.1579E-04, 4.5139E-07, 200, 1500, 1,2, 0.0255, 0.0883, 1.5635],
    [691, 'C4H10O', 'sec-butanol', 78_92_2, -0.07557, 4.1837E-04, 4.5288E-07, 200, 1500, 1,2, 0.0263, 0.0894, 1.5710],
    [692, 'C4H10O', 'tert-butanol', 75_65_0, -0.07652, 4.1543E-04, 4.5316E-07, 200, 1500, 1,2, 0.0248, 0.0876, 1.5662],
    [1215, 'C5H12O', '1-pentanol', 71_41_0, -0.07320, 3.6231E-04, 4.1679E-07, 200, 1500, 1,2, 0.0160, 0.0719, 1.4080],
    [1216, 'C5H12O', '2-pentanol', 6032_29_7, -0.07323, 3.6793E-04, 4.1180E-07, 200, 1500, 1,2, 0.0169, 0.0731, 1.4052],
    [1217, 'C5H12O', '3-pentanol', 584_02_1 , -0.07434, 3.7186E-04, 4.1469E-07, 200, 1500, 2, 0.0167, 0.0734, 1.4165],
    [1218, 'C5H12O', '2-methyl-1-butanol', 137_32_6, -0.07340, 3.6622E-04, 4.1865E-07, 200, 1500, 2, 0.0167, 0.0730, 1.4179],
    [1235, 'C5H12O', 'pentanol', 30899_19_5, -0.05353, 3.1574E-04, 4.2654E-07, 200, 1500, 2, 0.0267, 0.0785, 1.3798]
]


yaws_alcohol_details = ['Number', 'Formula', 'Name', 'CAS_Number', 'A', 'B', 'C', 'TMIN', 'TMAX', 'code', 'D_TMIN', 'D_25C', 'D_TMAX']
yaws_alcohols = {}

for i, alcohol in enumerate(yaws_alcohol_properties):
    yaws_alcohols[alcohol[2]] = dict(zip(yaws_alcohol_details, yaws_alcohol_properties[i]))

colours_list = ['r', 'k', 'b', 'limegreen', 'magenta']
yaws_alcohol_diffusion_fits = {}
for name in ['methyl_alcohol', 'ethyl_alcohol', 'propyl_alcohol', 'butanol', 'pentanol']:
    f = vapour_binary_diffusion_coefficient_func
    T = np.linspace(T_freezing - 10, T_freezing + 100 , 1000)
    D = yaws_diffusion_polynomial(T, yaws_alcohols[name]['A'],
                                  yaws_alcohols[name]['B'],
                                  yaws_alcohols[name]['C'])

    yaws_alcohol_diffusion_fits[name], _ = curve_fit(f, T, D)
# -

# Sanity check the diffusion fits by plotting them below:

if __name__ == '__main__':
    T_C = np.linspace(0, 100, 101)
    T_K = T_C + T_freezing

    plt.figure()
    for name in yaws_alcohols:
        D = yaws_diffusion_polynomial(T_K,
                                      yaws_alcohols[name]['A'],
                                      yaws_alcohols[name]['B'],
                                      yaws_alcohols[name]['C'])

        plt.plot(T_C, D, label=name)
        plt.legend(ncol=2)
        plt.title('All alcohols')
        plt.xlabel('T (℃)')
        plt.ylabel('$D_\infty$ (m$^2$/s)')
        plt.title('All alcohols, Yaws parameterisations')

    plt.figure()
    for name, colour in zip(['methyl_alcohol', 'ethyl_alcohol', 'propyl_alcohol', 'butanol', 'pentanol' ], colours_list):
        D = vapour_binary_diffusion_coefficient_func(T_K, *yaws_alcohol_diffusion_fits[name])
        plt.plot(T_C, D, '--', c=colour, label=(str(name) + ' fit'))

    plt.legend(ncol=2)
    plt.xlabel('T (℃)')
    plt.ylabel('$D_\infty$ (m$^2$/s)')
    plt.title('Selected alcohols with fits')

    plt.figure()
    for name, colour in zip(['methyl_alcohol', 'ethyl_alcohol', 'propyl_alcohol', 'butanol', 'pentanol' ], colours_list):
        D1 = vapour_binary_diffusion_coefficient_func(T_K, *yaws_alcohol_diffusion_fits[name])
        D2 = yaws_diffusion_polynomial(T_K,
                                       yaws_alcohols[name]['A'],
                                       yaws_alcohols[name]['B'],
                                       yaws_alcohols[name]['C'])
        plt.plot(T_C, D1 - D2, '-.', c=colour, label=(str(name) + ' fit - reference'))

    plt.legend()
    plt.xlabel('T (℃)')
    plt.ylabel('$D_\infty$ (m$^2$/s)')
    plt.title('Selected alcohols fit residuals');

# #### 2.1.4.3. Properties of pure ethannol

# +
molar_mass_ethanol = 46.07 # g/mol

def density_ethanol(T):
    """Density of ethanol as function of temperature.
   
    Functional form and fit parameters taken from (T_min, T_max = 191, 513):
    http://ddbonline.ddbst.de/DIPPR105DensityCalculation/DIPPR105CalculationCGI.exe?component=2-ethanol

    Args:
        T: temperature in Kelvin.
    Returns:
        The ethanol density in kg/m^3.
    """

    A, B, C, D = 99.3974, 0.310729, 513.18, 0.30514
    return rackett_equation(A, B, C, D, T)

def specific_heat_capacity_ethanol(T):
    """Specific heat capacity of ethanol vs temperature.

    Lacking better data, we assume this is constant wrt temperature taking
    numerical values from NIST:
    https://webbook.nist.gov/cgi/cbook.cgi?ID=C64175&Mask=2

    Args:
        T: temperature in Kelvin (though not currently used).
    Returns:
        Specific heat capacity in J/K/kg.
    """

    molar_heat_capacity_ethanol = 112.3 # J/K/mol
    return molar_heat_capacity_ethanol / molar_mass_ethanol * 1000 # J/K/kg

def specific_latent_heat_ethanol(T):
    """Enthalpy of vaporisation per unit mass.

    Functional form and parameters taken from Majer and Svoboda (1985), as referenced by NIST here:
    https://webbook.nist.gov/cgi/cbook.cgi?ID=C64175&Mask=4.

    Fit parameters for temperature range 298-469 K:
        A (kJ/mol)  50.43
        alpha       -0.4475
        beta         0.4989
        Tc (K)     513.9

    Returns:
        Specific latent heat in J/kg.
    """

    A, alpha, beta, T_c = 50.43, -0.4475, 0.4989, 513.9
    H_vap = enthalpy_vapourisation(50.43, -0.4475, 0.4989, 513.9, T) # J / mol
    return 1e-3 * H_vap / molar_mass_ethanol # J / kg

def equilibrium_vapour_pressure_ethanol(T):
    """Parameterisation of vapour pressure for ethanol vs temperature.

    Parameterisation taken from:
    https://webbook.nist.gov/cgi/inchi?ID=C71238&Mask=4&Type=ANTOINE&Plot=on

    This particular fit is probably a bit inaccurate below 20 deg c.

    Args:
       T: temperature in Kelvin.
    Returns:
       Vapour pressure in Pa.
    """

    P_bar = antoine_equation(T, 5.37229, 1670.409, -40.191)
    # Parameterisation above give pressure in bars, so we convert into Pa:
    return 1e5 * P_bar

def surface_tension_ethanol(T):
    """Surface tension of ethanol as a function of temperature.

    Functional form and parameters taken from:
    http://ddbonline.ddbst.de/DIPPR106SFTCalculation/DIPPR106SFTCalculationCGI.exe

    DIPPR106 parameters:
        A     131.38
        B       5.5437
        C      -8.4826
        D       4.3164
        E       0
        Tc    516.2
        Tmin  180
        Tmax  513

    Args:
        T: temperature in Kelvin.
    Returns:
        Surface tension in N/m.
    """

    return surface_tension(131.38, 5.5437, -8.4826, 4.3164, 0, 516.2, T)

Ethanol = Solvent(molar_mass_ethanol,
                  density_ethanol,
                  specific_heat_capacity_ethanol,
                  specific_latent_heat_ethanol,
                  equilibrium_vapour_pressure_ethanol,
                  VapourBinaryDiffusionCoefficient(0.2190e-4, T_freezing, 1.81),
                  surface_tension_ethanol)
# -

# Sanity check the ethanol properties by plotting them below:

if __name__ == '__main__':
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

    plt.figure()
    plt.plot(T_C, Ethanol.surface_tension(T_C + T_freezing), c = 'k')
    plt.xlabel('T (℃)')
    plt.ylabel('Surface tension / N / m');

# #### 2.1.4.4. Properties of pure propanol

# +
molar_mass_propanol = 60.0952 # g/mol

def density_propanol(T):
    """Density of propanol as function of temperature.
   
    Functional form and fit parameters taken from (T_min, T_max = 186, 507):
    http://ddbonline.ddbst.de/DIPPR105DensityCalculation/DIPPR105CalculationCGI.exe?component=2-propanol

    Args:
        T: temperature in Kelvin.
    Returns:
        The propanol density in kg/m^3.
    """

    A, B, C, D = 74.5237, 0.27342, 508.3, 0.235299
    return rackett_equation(A, B, C, D, T)

def specific_heat_capacity_propanol(T):
    """Specific heat capacity of propanol vs temperature.

    Lacking better data, we assume this is constant wrt temperature taking
    numerical values from NIST:
    https://webbook.nist.gov/cgi/cbook.cgi?ID=C71238&Mask=2

    Args:
        T: temperature in Kelvin (though not currently used).
    Returns:
        Specific heat capacity in J/K/kg.
    """

    molar_heat_capacity_propanol = 145 # J/K/mol
    return molar_heat_capacity_propanol / molar_mass_propanol * 1000 # J/K/kg

def specific_latent_heat_propanol(T):
    """Enthalpy of vaporisation per unit mass.

    Functional form and parameters taken from Majer and Svoboda (1985), as referenced by NIST here:
    https://webbook.nist.gov/cgi/cbook.cgi?ID=C64175&Mask=4.

    Fit parameters for temperature range 298-390 K:
        A (kJ/mol)  52.06
        alpha       -0.8386
        beta         0.6888
        Tc (K)     536.7

    Returns:
        Specific latent heat in J/kg.
    """

    A, alpha, beta, T_c = 52.06, -0.8386, 0.6888, 536.7
    H_vap = enthalpy_vapourisation(50.43, -0.4475, 0.4989, 513.9, T) # J / mol
    return 1e-3 * H_vap / molar_mass_propanol # J / kg

def equilibrium_vapour_pressure_propanol(T):
    """Parameterisation of vapour pressure for propanol vs temperature.

    Parameterisation taken from:
    https://webbook.nist.gov/cgi/inchi?ID=C71238&Mask=4&Type=ANTOINE&Plot=on

    This particular fit is probably a bit inaccurate below 20 deg c.

    Args:
       T: temperature in Kelvin.
    Returns:
       Vapour pressure in Pa.
    """

    P_bar = antoine_equation(T, 5.31, 1690.86, -51.80)
    # Parameterisation above give pressure in bars, so we convert into Pa:
    return 1e5 * P_bar

def surface_tension_propanol(T):
    """Surface tension of propanol as a function of temperature.

    Functional form and parameters taken from:
    http://ddbonline.ddbst.de/DIPPR106SFTCalculation/DIPPR106SFTCalculationCGI.exe

    DIPPR106 parameters:
        A      46.507
        B       0.90053
        C       0
        D       0
        E       0
        Tc    508.3
        Tmin  287
        Tmax  353

    Args:
        T: temperature in Kelvin.
    Returns:
        Surface tension in N/m.
    """

    return surface_tension(46.507, 0.90053, 0, 0, 0, 508.3, T)

Propanol = Solvent(molar_mass_propanol,
                density_propanol,
                specific_heat_capacity_propanol,
                specific_latent_heat_propanol,
                equilibrium_vapour_pressure_propanol,
                VapourBinaryDiffusionCoefficient(0.2190e-4, T_freezing, 1.81),
                surface_tension_propanol)
# -

# Sanity check the propanol properties by plotting them below:

if __name__ == '__main__':
    T_C = np.linspace(0, 100, 101)
    T_K = T_C + T_freezing
    fig, ((ax1,ax2),(ax3,a4)) = plt.subplots(ncols=2, nrows=2)

    ax1.plot(T_C, Propanol.density(T_K))
    ax1.set_xlabel('T (℃)')
    ax1.set_ylabel('density (kg/m3)')

    ax2.plot(T_C, Propanol.specific_latent_heat_vaporisation(T_K))
    #ax2.axhline(y=2264.705, ls='dashed')
    ax2.set_xlabel('T (℃)')
    ax2.set_ylabel('L (kJ/kg)')

    ax3.plot(T_C, Propanol.equilibrium_vapour_pressure(T_K))
    ax3.set_xlabel('T (℃)')
    ax3.set_ylabel('equilibrium vapour pressure (Pa)')

    ax4.plot(T_C, Propanol.vapour_binary_diffusion_coefficient(T_K))
    ax4.set_xlabel('T (℃)')
    ax4.set_ylabel('$D_\infty$ (m$^2$/s)')

    plt.figure()
    plt.plot(T_C, Propanol.surface_tension(T_K))
    plt.xlabel('T (℃)')
    plt.ylabel('Surface tension / N / m');

# #### 2.1.4.5. Properties of pure butanol

# +
molar_mass_butanol = 74.12 # g/mol

def density_butanol(T):
    """Density of butanol as function of temperature.

    Functional form and fit parameters taken from (T_min, T_max = 213, 558):
    http://ddbonline.ddbst.de/DIPPR105DensityCalculation/DIPPR105CalculationCGI.exe

    Args:
        T: temperature in Kelvin.
    Returns:
        The butanol density in kg/m^3.
    """

    A, B, C, D =  9.87035, 0.0998069, 568.017, 0.126276
    return rackett_equation(A, B, C, D, T)

def specific_heat_capacity_butanol(T):
    """Specific heat capacity of butanol vs temperature.

    Lacking better data, we assume this is constant wrt temperature taking
    numerical values from NIST:
    https://webbook.nist.gov/cgi/cbook.cgi?ID=C71238&Mask=2

    Args:
        T: temperature in Kelvin (though not currently used).
    Returns:
        Specific heat capacity in J/K/kg.
    """

    molar_heat_capacity_butanol = 177 # J/K/mol
    return molar_heat_capacity_butanol / molar_mass_butanol * 1000 # J/K/kg

def specific_latent_heat_butanol(T):
    """Enthalpy of vaporisation per unit mass.

    Functional form and parameters taken from Majer and Svoboda (1985), as referenced by NIST here:
    https://webbook.nist.gov/cgi/cbook.cgi?ID=C78922&Mask=4

    Fit parameters for temperature range 298 - 372 K:
        A (kJ/mol)  52.6
        alpha       -1.462
        beta         1.0701
        Tc (K)     536.

    Args:
        T: temperature in Kelvin.
    Returns:
        Specific latent heat in J/kg.
    """

    A, alpha, beta, T_c = 52.6, -1.462, 1.0701, 536
    H_vap = enthalpy_vapourisation(50.43, -0.4475, 0.4989, 513.9, T) # J / mol
    return 1e-3 * H_vap / molar_mass_butanol # J / kg

def equilibrium_vapour_pressure_butanol(T):
    """Parameterisation of vapour pressure for butanol vs temperature.

    Parameterisation taken from:
    https://webbook.nist.gov/cgi/cbook.cgi?ID=C71363&Mask=4&Type=ANTOINE&Plot=on

    This particular fit is probabaly a bit inaccurate below 20 deg c.

    Args:
       T: temperature in Kelvin.
    Returns:
       Vapour pressure in Pa.
    """

    P_bar = antoine_equation(T, 4.55, 1351.56, -93.34)
    # Parameterisation above give pressure in bars, so we convert into Pa:
    return 1e5 * P_bar

def surface_tension_butanol(T):
    """Surface tension of ethanol as a function of temperature.

    Functional form and parameters taken from:
    http://ddbonline.ddbst.de/DIPPR106SFTCalculation/DIPPR106SFTCalculationCGI.exe

    DIPPR106 parameters:
        A     72.697
        B       3.0297
        C      -4.2681
        D       2.4776
        E       0
        Tc    562.9
        Tmin  238
        Tmax  543

    Args:
        T: temperature in Kelvin.
    Returns:
        Surface tension in N/m.
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

# Sanity check the butanol properties by plotting them below:

if __name__ == '__main__':
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

    plt.figure()
    plt.plot(T_C, Butanol.surface_tension(T_C + T_freezing), c = 'g')
    plt.xlabel('T (℃)')
    plt.ylabel('Surface tension / N / m')
    plt.show()
