# Useful functions for extracting E-AIM data

def rename_eaim_columns(df):
    """Takes a df with the headings from E-AIM and makes the columns nice for using in pandas.

    Gets rid of (aq), (g) and spaces.

    Args:
        df: the table (assumed pandas dataframe).
    Returns:
        A list of columns.
    """
    columns = df.columns.str.strip()
    columns = columns.str.replace(' ', '_')
    columns = columns.str.replace('(aq)', '', regex=False)
    columns = columns.str.replace('(g)', '', regex=False)
    return columns

def mfs_from_molality (molality_series, solute_molar_mass):
    """Takes molality data of a species - *THIS MUST MATCH THE MOLALITY OF THE SOLUTE*
        i.e. n_species = n_solute. e.g. for MgCl_2, must use molality_Mg.

    Takes molar mass of solute in g/mol.

    Args:
        molality_series: ?
        solute_molar_mass: units?
    Returns:
        Mass fraction series.
    """
    return 1 / ( 1 + (1 / (molality_series * solute_molar_mass * 1e-3)))

def eaim_water_activity(df):
    """Water activity from reference data (from E-AIM).

    Args:
        df: the table of E-AIM data (assumed pandas dataframe).
    """
    return df.f_H2O * df.x_H2O

def fit_eaim_activity(mfs, activity, degree):
    """Perform a fit of solvent activity (a) vs mass fraction of solute (mfs) on some reference data.

    We constrain the fit so that a(mfs=0) = 1 and a(mfs=1) = 0, the physically correct limits.
    The resulting polynomial fit then has the form:

        a = c_0 + c_1 * mfs + c_2 * mfs**2 + ... + c_(n-1) * mfs**(n-1) + c_n * mfs**n

    where c_0 = 1 and c_n = -\sum_{i=0}^{n-1} c_i in order to have the correct limits.

    Args:
        mfs: reference mass fraction of solute.
        activity: reference solvent activity.
        degree: degree of polynomial fit.
    Returns:
        The fit function.
    """

    coefficients = lambda x: np.concatenate([[-1-np.sum(x)], x, [1]])
    fit_func = lambda mfs,*x: np.polyval(coefficients(x), mfs)
    fit = curve_fit(fit_func, mfs, activity, p0=np.zeros(degree-1))[0]
    return lambda mfs: fit_func(mfs, *fit)

import os

location = os.path.dirname(os.path.realpath(__file__))

import pandas
aqueous_NaCl_data_path = os.path.join(location, 'EDB data for benchmarking', 'NaCl_all.csv')
aqueous_NaCl_data = pandas.read_csv(aqueous_NaCl_data_path)
aqueous_NaCl_data.columns = rename_eaim_columns(aqueous_NaCl_data)

from chemicals import periodic_table
NaCl_molar_mass = periodic_table.Na.MW + periodic_table.Cl.MW

aqueous_NaCl_data['mass_fraction_solute'] = mfs_from_molality(aqueous_NaCl_data.m_Na, NaCl_molar_mass)
aqueous_NaCl_data['solvent_activity'] = eaim_water_activity(aqueous_NaCl_data)
