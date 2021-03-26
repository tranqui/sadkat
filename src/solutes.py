# + ignore="True"
from solvents import *
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
if __name__ == '__main__':
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
