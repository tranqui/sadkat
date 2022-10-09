# + ignore="True"
from sadkat.solvents import *
# -

# ## 2.2. Solutes
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
    x = np.real(np.min(x[(np.isreal(x)) & (x >= 0)]))
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

if __name__ == '__main__':
    # Dummy parameters for the fits.
    forward_fit = ActivityVsMfsParameterisation([1, -1, 1, -1])
    backward_fit = MfsVsActivityParameterisation([1, -1, 1, -1])

    mfs = np.linspace(0, 1, 100)
    aw = np.linspace(0, 1, 100)

    fig, (ax1, ax2) = plt.subplots(ncols=2)

    ax1.plot(mfs, forward_fit.solvent_activity_from_mass_fraction_solute(mfs), label='regular fit')
    ax1.plot(forward_fit.mass_fraction_solute_from_solvent_activity(aw), aw, ls='--', label='inverted fit')
    ax1.legend(loc='best')
    ax1.set_title('forward parameterisation')

    ax2.plot(backward_fit.mass_fraction_solute_from_solvent_activity(aw), aw, ls='--', label='regular fit')
    ax2.plot(mfs, backward_fit.solvent_activity_from_mass_fraction_solute(mfs), label='inverted fit')
    ax2.legend(loc='best')
    ax2.set_title('backward parameterisation')

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
    solvent_activity_fit: object=None  # a fitting function for solvent activity (unitless) in
                                       # terms of the mass fraction of solute, or None (the
                                       # default) to assume Raoult's law for ideal mixtures

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
        """Invert solvent activity to find the mass fraction of the solute at a particular
        solvent activity.

        This is useful for determining the equilibrium state of a droplet by matching
        vapour pressures at the liquid-gas boundary.

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
            The density at this value of mfs (kg/m^3).
        """
        return 1 / (mfs / self.pure_solute_density + (1-mfs) / self.pure_solvent_density)
# -

# Parameterisations for some common solutions:

# +
from chemicals import periodic_table
NaCl_molar_mass = periodic_table.Na.MW + periodic_table.Cl.MW
aqueous_NaCl = Solution(Water, 58.44, 2, 0.3, 2170,
                        DensityVsMassFractionFit([998.2 , -55.33776, 1326.69542, -2131.05669, 2895.88613, -940.62808]),
                        ActivityVsMfsParameterisation(np.flipud([48.5226539, -158.04388699, 186.59427048, -93.88696437, 19.28939256, -2.99894206, -0.47652352, 1.])))

aqueous_NaCl_a_ideal = Solution(Water, 58.44, 2, 0.3, 2170,
                                DensityVsMassFractionFit([998.2 , -55.33776, 1326.69542, -2131.05669, 2895.88613, -940.62808]),
                                ActivityVsMfsParameterisation(np.flipud([48.5226539, -158.04388699, 186.59427048, -93.88696437, 19.28939256, -2.99894206, -0.47652352, 1.])))

aqueous_NaCl_d_linear = Solution(Water, 58.44, 2, 0.3, 2170,
                                 DensityVsMassFractionFit([998.2 ,1162]),
                                 ActivityVsMfsParameterisation(np.flipud([48.5226539, -158.04388699, 186.59427048, -93.88696437, 19.28939256, -2.99894206, -0.47652352, 1.])))

aqueous_NaCl_d_half = Solution(Water, 58.44, 2, 0.3, 2170,
                               DensityVsMassFractionFit([998.2 , 0.5 * -55.33776, 0.5 * 1326.69542, 0.5 * -2131.05669, 0.5 * 2895.88613, 0.5 * -940.62808]),
                               ActivityVsMfsParameterisation(np.flipud([48.5226539, -158.04388699, 186.59427048, -93.88696437, 19.28939256, -2.99894206, -0.47652352, 1.])))

# Other inorganic salts

aqueous_KCl = Solution(Water, 74.5513, 2, 0.25, 1980,
                       DensityVsMassFractionFit([996.75598478, 28.58359584, 471.74984825, 324.24333223]),
                       ActivityVsMfsParameterisation(np.flipud( [ 9.38399042e+03, -5.29769730e+04, 1.17504454e+05, -1.21998524e+05, 4.09435471e+04, 3.25878163e+04, -3.45994305e+04, 5.81753302e+03, 6.89074445e+03, -4.69316888e+03, 1.31490992e+03, -1.87352995e+02, 1.21808905e+01, -7.26450632e-01, 1.00000000e+00])))

aqueous_NaI = Solution(Water, 149.89, 2, 0.65, 3670,
                       DensityVsMassFractionFit([992.70440241, 119.97468447, 308.3296718 , 977.06460272]),
                       ActivityVsMfsParameterisation(np.flipud([-4.58977533e+03, 2.36936878e+04, -5.04825135e+04, 5.64293503e+04, -3.40855148e+04, 9.50188876e+03, -2.49232286e+02, 1.20624905e+02, -5.53378784e+02, 2.59324280e+02, -4.65900167e+01, 1.36924533e+00, 6.38800214e-04, -2.41266758e-01, 1.00000000e+00])))

aqueous_KI = Solution(Water, 166.0028, 2, 0.6, 3120,
                      DensityVsMassFractionFit([988.83030777,  158.42666625,  -32.49513622, 1256.6891563]),
                      ActivityVsMfsParameterisation(np.flipud([-4.95757159e+03, 1.99508232e+04, -2.87027928e+04, 1.23204357e+04, 1.07871868e+04, -1.39672501e+04, 4.01678577e+03, 1.60042543e+03, -1.31354702e+03, 2.61615629e+02, 1.39094494e+01, -1.19736506e+01, 1.19828942e+00, -2.44979095e-01, 1.00000000e+00])))

aqueous_LiNO3 = Solution(Water, 68.946, 2, 0.34, 2380,
                         DensityVsMassFractionFit([998.56832988, -60.79985016, 1005.61549769, -1085.10992993, 1177.67900264]),
                         ActivityVsMfsParameterisation(np.flipud([-4.23252154e+03, 2.17183817e+04, -4.34823030e+04, 3.84652543e+04, -4.46004526e+03, -1.98879519e+04, 1.58595251e+04, -2.97181926e+03, -2.16742018e+03, 1.52132958e+03, -4.12257237e+02, 5.37352358e+01, -4.51076657e+00, -3.96785308e-01, 1.00000000e+00])))

aqueous_NaNO3 = Solution(Water, 84.9947, 2, 0.48, 2260,
                         DensityVsMassFractionFit([990.53882439, 181.01413995, -494.89819949, 2895.56981334, -2999.14475279, 1526.74834627]),
                         ActivityVsMfsParameterisation(np.flipud([-3.59764720e+03, 1.72951494e+04, -3.34613879e+04, 3.18176263e+04, -1.33409555e+04, 4.79558719e+02, -6.59464410e+02, 3.23488952e+03, -2.47664569e+03, 8.20527382e+02, -1.11013532e+02, -3.11492521e+00, 1.99338440e+00, -5.15557822e-01, 1.00000000e+00])))

aqueous_KNO3 = Solution(Water, 101.1032, 2, 0.28, 2110,
                        DensityVsMassFractionFit([997.1673202, 23.13666049, 473.87886124, 339.61281717]),
                        ActivityVsMfsParameterisation(np.flipud([-8.01841811e+02, 3.13053769e+03, -4.26110013e+03, 1.64775304e+03, 1.22569565e+03, -9.07034756e+02, -4.36504870e+02, 4.65383398e+02, 4.83605721e+01, -1.76260611e+02, 7.87046110e+01, -1.60383987e+01, 1.70555883e+00, -3.59952132e-01,  1.00000000e+00])))

aqueous_RbNO3 = Solution(Water, 147.473, 2, 0.39, 3110,
                         DensityVsMassFractionFit([995.49397656, 85.49864715, 81.33053489, 991.79067905]),
                         ActivityVsMfsParameterisation(np.flipud([-4.39043562, 5.41778218, -2.08008163,  0.23378074, -0.18104567, 1.])))

aqueous_Na2SO4 = Solution(Water, 142.04, 3, 0.3, 2660,
                          DensityVsMassFractionFit([9.98003117e+02, -8.65783619e+01, 4.30749205e+03, -6.52155268e+04, 4.97341694e+05, -2.05621390e+06, 5.05960147e+06, -7.64562419e+06, 6.97822368e+06, -3.53493475e+06, 7.64101736e+05]),
                          ActivityVsMfsParameterisation(np.flipud([3.10491640e+03, -1.44633278e+04, 2.68928321e+04, -2.45145530e+04, 9.77420691e+03, 6.35871766e+02, -2.05298685e+03, 7.51626505e+02, -1.46748722e+02, -2.89668158e+00, 3.55591961e+01, -1.87591912e+01, 3.63716532e+00, -3.77883865e-01, 1.00000000e+00,])))

aqueous_K2SO4 = Solution(Water, 174.259, 3, 0.11, 2660,
                         DensityVsMassFractionFit([997.84866879, 8.23246363, 734.54615575, 210.64816357]),
                         ActivityVsMfsParameterisation(np.flipud([-1.88138199e+04, 1.03068495e+05, -2.30116240e+05, 2.59420208e+05, -1.35783431e+05, -1.76822793e+03, 3.89160605e+04, -1.51183352e+04, -2.22686968e+03, 3.33557836e+03, -1.05936351e+03, 1.55377081e+02, -1.04407278e+01, 9.55263257e-03,  1.00000000e+00])))

aqueous_MgSO4 = Solution(Water, 120.366, 2, 0.17, 2660,
                         DensityVsMassFractionFit([995.05558478, 63.40780942, 634.1079356,  780.69806718]),
                         ActivityVsMfsParameterisation(np.flipud([-4.59916396e+02, 2.85037407e+03, -6.28221570e+03, 5.52474650e+03, -1.06307021e+02, -2.80056023e+03, 8.45576689e+02, 1.19994793e+03, -1.10241166e+03, 3.97510492e+02, -6.70379431e+01, -1.16460310e+00, 6.20603120e-01, -1.62727222e-01, 1.00000000e+00])))

## "real world" solutions
# Warning: besides the density fits, the parameters below for DLF/AS are just
# dummy placeholder values as far as I am aware.

aqueous_DLF  = Solution(Water,  1.00, 2, 1.0, 1700,
                        DensityVsMassFractionFit([995.78, 262.92   , -606.15   ,  1135.53   ,    0      ,    0]))

aqueous_AS   = Solution(Water,  1.00, 2, 1.0, 1200,
                        DensityVsMassFractionFit([997   , -43.88148,  397.75347,  -100.99474,    0      ,    0]))

# Solutions using different solvents

#dummy solute properties - USE ONLY AS PURE SOLVENT - mfs = 0
ethanolic_NaCl = Solution(Ethanol, 58.44, 2, 0.3, 2170,
                          DensityVsMassFractionFit([789]))

#dummy solute properties - USE ONLY AS PURE SOLVENT
propanolic_NaCl = Solution(Propanol, 58.44, 2, 0.3, 2170,
                          DensityVsMassFractionFit([803]))

#dummy solute properties - USE ONLY AS PURE SOLVENT
butanolic_NaCl = Solution(Butanol, 58.44, 2, 0.3, 2170,
                          DensityVsMassFractionFit([810]))

all_solutions = {'NaCl in water': aqueous_NaCl,
                 'KCl in water': aqueous_KCl,
                 'NaI in water': aqueous_NaI,
                 'KI in water': aqueous_KI,
                 'LiNO3 in water': aqueous_LiNO3,
                 'NaNO3 in water': aqueous_NaNO3,
                 'KNO3 in water': aqueous_KNO3,
                 'RbNO3 in water': aqueous_RbNO3,
                 'Na2SO4 in water': aqueous_Na2SO4,
                 'K2SO4 in water': aqueous_K2SO4,
                 'MgSO4 in water': aqueous_MgSO4}

                 # 'DLF in water': aqueous_DLF,
                 # 'AS in water': aqueous_AS,
                 # 'NaCl in ethanol': ethanolic_NaCl,
                 # 'NaCl in propanol': propanolic_NaCl,
                 # 'NaCl in butanol': butanolic_NaCl}
# -

# Identical versions of the above solutions but with a simpler density parameterisation (volume additivity):

from copy import copy
all_solutions_volume_additivity = {}
for label, solution in all_solutions.items():
    additive_solution = copy(solution)
    rho1, rho2 = solution.density(0), solution.solid_density
    additive_solution.density = VolumeAdditivityFit(rho1, rho2)
    all_solutions_volume_additivity[label] = additive_solution

# Test parameterisations of aqueous NaCl solution against verified data:

if __name__ == '__main__':

    from sadkat.eaim import aqueous_NaCl_data as NaCl_data_EAIM

    mfs = np.linspace(0, 1, 100)
    fig, (ax1, ax2) = plt.subplots(ncols=2)

    label = 'NaCl in water'
    additive_solute = all_solutions_volume_additivity[label]

    pl1, = ax1.plot(mfs, aqueous_NaCl.density(mfs), zorder=0, label='Fit')
    pl2, = ax1.plot(NaCl_data_EAIM.mass_fraction_solute, 1000 * NaCl_data_EAIM.Density, '.', label='E-AIM')

    ax1.set_xlabel('MFS')
    ax1.set_ylabel('density (kg/m$^3$)')
    ax1.legend()

    aw = np.linspace(0, 1, 100)

    ax2.plot(mfs, aqueous_NaCl.solvent_activity(mfs), c=pl1.get_color(), zorder=0, label='Fit')
    ax2.plot(NaCl_data_EAIM.mass_fraction_solute, NaCl_data_EAIM.solvent_activity, '.', c=pl2.get_color(), label='E-AIM')

    ax2.set_xlabel('MFS')
    ax2.set_ylabel('solvent activity')
    ax2.legend()

# Sanity check the parameterisations of the solutions by plotting some of their properties below:

if __name__ == '__main__':
    mfs = np.linspace(0, 1, 100)
    fig, (ax1, ax2) = plt.subplots(ncols=2)

    for label, solute in all_solutions.items():
        additive_solute = all_solutions_volume_additivity[label]

        pl, = ax1.plot(mfs, solute.density(mfs), label=label)
        ax1.plot(mfs, additive_solute.density(mfs), '--', c=pl.get_color(), label=('%s (volume additivity)' % label))
        ax1.plot(1, solute.solid_density, 'o', mfc='None', c=pl.get_color())

    ax1.set_xlabel('MFS')
    ax1.set_ylabel('density (kg/m$^3$)')
    ax1.legend(loc='upper center', bbox_to_anchor = (0.5,-0.2))

    aw = np.linspace(0, 1, 100)
    for label, solute in all_solutions.items():
        pl, = ax2.plot(mfs, solute.solvent_activity(mfs), label=label)
        ax2.plot(solute.mass_fraction_solute_from_solvent_activity(aw), aw, '--', c=pl.get_color(), label=('%s (inverted fit)' % label))

    ax2.set_xlabel('MFS')
    ax2.set_ylabel('solvent activity')
    ax2.legend(loc='upper center', bbox_to_anchor = (0.5,-0.2))

# Test various ways to parameterise aqueous NaCl data:

if __name__ == '__main__':

    mfs = np.linspace(0, 1, 100)
    fig, (ax1, ax2) = plt.subplots(ncols=2)

    label = 'NaCl in water'
    additive_solute = all_solutions_volume_additivity[label]

    pl1, = ax1.plot(mfs, aqueous_NaCl.density(mfs), zorder=0, label='E-AIM Fit')
    pl2, = ax1.plot(NaCl_data_EAIM.mass_fraction_solute, 1000 * NaCl_data_EAIM.Density, '.', label='E-AIM Data')

    pl3, = ax1.plot(mfs, aqueous_NaCl_d_linear.density(mfs), ls='--', zorder=0, label='Linear in $\sqrt{MFS}$')
    pl4, = ax1.plot(mfs, aqueous_NaCl_d_half.density(mfs), ls='-.', zorder=0, label=r'$\frac{1}{2}$ $\times$ E-AIM Fit')

    ax1.set_xlabel('MFS')
    ax1.set_ylabel('density (kg/m$^3$)')
    ax1.legend()

    aw = np.linspace(0, 1, 100)

    ax2.plot(mfs, aqueous_NaCl.solvent_activity(mfs), c=pl1.get_color(), zorder=0, label='E-AIM Fit')
    ax2.plot(NaCl_data_EAIM.mass_fraction_solute, NaCl_data_EAIM.solvent_activity, '.', c=pl2.get_color(), label='E-AIM Data')

    ax2.plot(mfs, aqueous_NaCl_a_ideal.solvent_activity(mfs), ls='--', c='r', zorder=0, label="Raoult's Law")

    ax2.set_xlabel('MFS')
    ax2.set_ylabel('solvent activity')
    ax2.legend()

    plt.show()
