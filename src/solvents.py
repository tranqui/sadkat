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
    """Takes fitting parameters and temperature (K) and returns surface tension (N/m).
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

if __name__ == '__main__':
    R_range = np.arange(1e-9,30e-9, 1e-10)

    curve_range = kelvin_effect(0.073, 997, 0.018, 293, 2300, R_range )

    plt.figure()
    plt.plot(R_range/1e-9,curve_range/2300)
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
    """
    T_C = T - T_freezing # Celsius
    return 1e3*0.61161 * np.exp((18.678 - (T_C / 234.5)) * (T_C / (257.14 + T_C))) # Pa


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

if __name__ == '__main__':
    def plot_solvent_properties(solvent, T = np.linspace(T_freezing, T_freezing + 99.9, 201)):
        fig, axes = plt.subplots(ncols = 2, nrows = 3, sharex = True, dpi = resolution, figsize = (figure_width, figure_height * 2))
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

        plt.show()

if __name__ == '__main__':
    plot_solvent_properties(Water)

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

if __name__ == '__main__':
    plt.figure()
    plt.plot(T_C, Water.surface_tension(T_C + T_freezing))
    plt.xlabel('T (℃)')
    plt.ylabel('Surface tension / N / m')

def rackett_equation(A, B, C, D, T):
    """ Rackett equation, takes coefficients and T (in Kelvin) and returns density.
    Taken from http://ddbonline.ddbst.de/DIPPR105DensityCalculation/DIPPR105CalculationCGI.exe?component=Methanol
    """

    density = (A) / (B ** (1 + ((1 - (T / C))) ** D))

    return density

def antoine_equation(T, A, B, C):
    """T in deg C, returns P in whatever parameters are for"""

    P =  10 ** ( A - (B / (T + C)))

    return P

# +
# NIST latent heat of vapourisation parameterisation
#
# ![image.png](attachment:image.png)
# -

def enthalpy_vapourisation(A, alpha, beta, T_c, T):
    """Rreceived T in K and coefficients. Returns H_vap in J/mol"""
    T_r = T / T_c

    H_vap = 1000 * ( A * np.exp( - alpha * T_r) * (1 - T_r) ** beta )

    return H_vap

def temp_K_to_C(temperature_K):
    return temperature_K - 273.15

# Alcohol diffusion coefficients
#
# Reference diffusion values: https://www.sciencedirect.com/science/article/pii/B978032328658900010X

# +
### Alcohol diffusion coefficients
#Reference Diffusion values https://www.sciencedirect.com/science/article/pii/B978032328658900010X
T_C = np.linspace(0, 100, 101)

def yaws_diffusion_polynomial(T,A,B,C):
    return 1e-4 * (A + B*T + C*T**2)

def vapour_binary_diffusion_coefficeint_func(T, D_ref, T_ref, lam):
    return D_ref * (T / T_ref)**lam

yaws_alcohol_props = [[40, 'CH4O', 'methyl_alcohol', 67_56_1, -0.12691, 7.2728E-04, 7.0516E-07, 200, 1500, 1,2, 0.0469, 0.1526, 2.5506],
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
                 [1235, 'C5H12O', 'pentanol', 30899_19_5, -0.05353, 3.1574E-04, 4.2654E-07, 200, 1500, 2, 0.0267, 0.0785, 1.3798]]


yaws_alcohol_details = ['Number', 'Formula', 'Name', 'CAS_Number', 'A', 'B', 'C', 'TMIN', 'TMAX', 'code', 'D_TMIN', 'D_25C', 'D_TMAX']
yaws_alcohols = {}

for i, alcohol in enumerate(yaws_alcohol_props):
    yaws_alcohols[alcohol[2]] = dict(zip(yaws_alcohol_details, yaws_alcohol_props[i]))

colours_list = ['r', 'k', 'b', 'limegreen', 'magenta']
yaws_alcohol_diffusion_fits = {}

for name in yaws_alcohols:
    plt.plot(T_C,
             yaws_diffusion_polynomial(T_C+T_freezing,
                     yaws_alcohols[name]['A'],
                     yaws_alcohols[name]['B'],
                     yaws_alcohols[name]['C']),
             label = name)
plt.legend(ncol = 2, fontsize = '14')
plt.title('All alcohols', fontsize = 24)
plt.xlabel('T (℃)')
plt.ylabel('$D_\infty$ (m$^2$/s)')
plt.title('All alcohols, Yaws parameterisations', fontsize = 24)
plt.show()

for name, colour in zip(['methyl_alcohol', 'ethyl_alcohol', 'propyl_alcohol', 'butanol', 'pentanol' ], colours_list):
    plt.plot(T_C,
             yaws_diffusion_polynomial(T_C+T_freezing,
                     yaws_alcohols[name]['A'],
                     yaws_alcohols[name]['B'],
                     yaws_alcohols[name]['C']),
             color = colour,
             label = name)

    yaws_alcohol_diffusion_fits[name], _ = curve_fit(vapour_binary_diffusion_coefficeint_func,
                                                     np.linspace(T_freezing - 10,T_freezing + 100 ,1000),
                                                     yaws_diffusion_polynomial(np.linspace(T_freezing - 10,T_freezing + 100 ,1000),
                                                                               yaws_alcohols[name]['A'],
                                                                               yaws_alcohols[name]['B'],
                                                                               yaws_alcohols[name]['C']))

for name, colour in zip(['methyl_alcohol', 'ethyl_alcohol', 'propyl_alcohol', 'butanol', 'pentanol' ], colours_list):
    plt.plot(T_C,
             vapour_binary_diffusion_coefficeint_func(T_C + T_freezing,
                                                      yaws_alcohol_diffusion_fits[name][0],
                                                      yaws_alcohol_diffusion_fits[name][1],
                                                      yaws_alcohol_diffusion_fits[name][2]),
             color = colour,
             linestyle = '--',
             label = str(name) + ' fit')



plt.legend(ncol = 2)
plt.xlabel('T (℃)')
plt.ylabel('$D_\infty$ (m$^2$/s)')
plt.title('Selected alcohols with fits', fontsize = 24)
plt.show()

for name, colour in zip(['methyl_alcohol', 'ethyl_alcohol', 'propyl_alcohol', 'butanol', 'pentanol' ], colours_list):
    plt.plot(T_C,
             vapour_binary_diffusion_coefficeint_func(T_C + T_freezing,
                                                      yaws_alcohol_diffusion_fits[name][0],
                                                      yaws_alcohol_diffusion_fits[name][1],
                                                      yaws_alcohol_diffusion_fits[name][2]) - yaws_diffusion_polynomial(T_C+T_freezing,
                                                                                                                        yaws_alcohols[name]['A'],
                                                                                                                        yaws_alcohols[name]['B'],
                                                                                                                        yaws_alcohols[name]['C']),
             color = colour,
             linestyle = '-.',
             label = str(name) + ' fit - reference' )
plt.legend()
plt.xlabel('T (℃)')
plt.ylabel('$D_\infty$ (m$^2$/s)')
plt.title('Selected alcohols fit residuals', fontsize = 24)

plt.show()

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

if __name__ == '__main__':
    plt.figure()
    plt.plot(T_C, Ethanol.surface_tension(T_C + T_freezing), c = 'k')
    plt.xlabel('T (℃)')
    plt.ylabel('Surface tension / N / m')

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

if __name__ == '__main__':
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

if __name__ == '__main__':
    plt.figure()
    plt.plot(T_C, Propanol.surface_tension(T_C + T_freezing), c = 'r')
    plt.xlabel('T (℃)')
    plt.ylabel('Surface tension / N / m')

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

if __name__ == '__main__':
    plt.figure()
    plt.plot(T_C, Butanol.surface_tension(T_C + T_freezing), c = 'g')
    plt.xlabel('T (℃)')
    plt.ylabel('Surface tension / N / m')
    plt.show()
