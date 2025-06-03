# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019--2025  Matthew P. Humphreys  (GNU GPLv3)
"""Calibrate and solve titration datasets.

This is the lowest level of the three layers of processing functions.
Titration data need to be imported and converted into separate NumPy arrays in
order to work with the functions here.

First-guess functions
---------------------
gran_function
gran_alkalinity
gran_emf0s
gran_guesses

Processing functions
--------------------
totals_ks
add_titrant_totals

Solver functions
----------------
solve_emf
solve_pH

Calibration functions
---------------------
calibrate_emf
calibrate_pH
"""

import warnings
from collections import namedtuple
from warnings import warn

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import linregress

from . import constants, convert, interface, simulate
from .settings import kwargs_least_squares


gar_vars = ("alkalinity", "gfunc", "used", "intercept_x", "lr")
GranAlkalinityResult = namedtuple("GranAlkalinityResult", gar_vars)
GranGuessesResult = namedtuple(
    "GranGuessesResult", (*gar_vars, "emf0", "emf0s", "pH")
)
SolveEmfResult = namedtuple(
    "SolveEmfResult",
    (
        "alkalinity",
        "emf0",
        "used",
        "pH",
        "alkalinity_std",
        "alkalinity_all",
        "opt_result",
        "ggr",
    ),
)
SolvePhResult = namedtuple(
    "SolvePhResult",
    (
        "alkalinity",
        "used",
        "alkalinity_std",
        "alkalinity_all",
    ),
)
# These are the kwargs that can be passed to `totals_ks`
keys_totals_ks = {
    "dic",
    "k_alpha",
    "k_ammonia",
    "k_beta",
    "k_bisulfate",
    "k_borate",
    "k_carbonic_1",
    "k_carbonic_2",
    "k_fluoride",
    "k_phosphoric_1",
    "k_phosphoric_2",
    "k_phosphoric_3",
    "k_silicate",
    "k_sulfide",
    "k_water",
    "opt_k_bisulfate",
    "opt_k_carbonic",
    "opt_k_fluoride",
    "opt_pH_scale",
    "opt_total_borate",
    "total_alpha",
    "total_ammonia",
    "total_beta",
    "total_borate",
    "total_fluoride",
    "total_phosphate",
    "total_silicate",
    "total_sulfate",
    "total_sulfide",
}


def gran_function(titrant_mass, emf, temperature, analyte_mass):
    """Calculate Gran-plot estimator (DAA03 eq. 10) using all provided data.

    Parameters
    ----------
    titrant_mass : array-like float
        Mass of titrant in kg.
    emf : array-like float
        EMF measured in titrant-analyte mixture in mV.
    temperature : array-like float
        Temperature of titrant-analyte mixture in °C.
    analyte_mass : float
        Mass of analyte in kg.

    Returns
    -------
    array-like float
        Gran-plot estimator (DAA03 eq. 10).
    """
    return (titrant_mass + analyte_mass) * np.exp(
        emf
        * constants.faraday
        / (constants.ideal_gas * (temperature + constants.absolute_zero))
    )


def gran_alkalinity(
    titrant_mass,
    emf,
    temperature,
    analyte_mass,
    titrant_molinity,
    titrant_normality=1,
):
    """Gran-plot first estimate of total alkalinity.

    Parameters
    ----------
    titrant_mass : array-like float
        Mass of titrant in kg.
    emf : array-like float
        EMF measured in titrant-analyte mixture in mV.
    temperature : array-like float
        Temperature of titrant-analyte mixture in °C.
    analyte_mass : float
        Mass of analyte in kg.
    titrant_molinity : float
        Molinity of titrant in mol/kg-sol.
    titrant_normality : float, optional
        Titrant normality, by default 1 (e.g., for HCl).

    Returns
    -------
    GranAlkalinityResult : namedtuple with the fields
        alkalinity : float
            Alkalinity estimate in mol/kg-sol.
        gfunc : array-like float
            Gran function values from `gran_function`.
        used : array-like bool
            Which Gran function points are used.
        intercept_x : float
            x-axis intercept computed from `lr`.
        lr : scipy.stats.LinregressResult
            Output from `scipy.stats.linregress`.
    """
    # Calculate Gran estimates
    gfunc = gran_function(titrant_mass, emf, temperature, analyte_mass)
    if (np.diff(gfunc) <= 0).any():
        warn("Gran function is not always increasing")
    # Determine which Gran estimates to use for fitting
    used = np.full(gfunc.shape, False)
    # use_logic = np.diff(gfunc) > np.max(np.diff(gfunc)) * 0.9
    use_logic = gfunc > 0.1 * np.max(gfunc)
    use_from = use_logic.nonzero()[0][0]
    used[use_from:] = True
    if used.sum() < 3:
        warn("Fewer than 3 data points available for linear regression")
    # Do linear regression
    lr = linregress(titrant_mass[used], gfunc[used])
    if lr.rvalue < 0.95:
        warn("Linear regression rvalue lower than 0.95")
    intercept_x = -lr.intercept / lr.slope
    alkalinity = (
        intercept_x * titrant_molinity * titrant_normality / analyte_mass
    )
    return GranAlkalinityResult(alkalinity, gfunc, used, intercept_x, lr)


def gran_emf0s(
    titrant_mass,
    emf,
    temperature,
    analyte_mass,
    titrant_molinity,
    alkalinity,
    titrant_normality=1,
):
    """Calculate Gran-plot estimator (DAA03 eq. 11) for EMF0.

    Parameters
    ----------
    titrant_mass : array-like float
        Mass of titrant in kg.
    emf : array-like float
        EMF measured in titrant-analyte mixture in mV.
    temperature : array-like float
        Temperature of titrant-analyte mixture in °C.
    analyte_mass : float
        Mass of analyte in kg.
    titrant_molinity : float
        Molinity of titrant in mol/kg-sol.
    alkalinity : float
        (Estimated) total alkalinty in mol/kg-sol.
    titrant_normality : float, optional
        Titrant normality, by default 1 (e.g., for HCl).

    Returns
    -------
    array-like float
        Gran EMF0 estimator (DAA03 eq. 11).
    """
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            category=RuntimeWarning,
            message="invalid value encountered in log",
        )
        emf0s = emf - (
            constants.ideal_gas
            * (temperature + constants.absolute_zero)
            / constants.faraday
        ) * np.log(
            (
                titrant_mass * titrant_molinity * titrant_normality
                - analyte_mass * alkalinity
            )
            / (titrant_mass + analyte_mass)
        )
    return emf0s


def gran_guesses(
    titrant_mass,
    emf,
    temperature,
    analyte_mass,
    titrant_molinity,
    titrant_normality=1,
):
    """Calculate Gran-plot first guesses for alkalinity, EMF0 and pH.

    Parameters
    ----------
    titrant_mass : array-like
        Mass of titrant in kg.
    emf : array-like
        EMF measured across the titrant-analyte mixture in mV.
    temperature : array-like
        Temperature of titrant-analyte mixture in °C.
    analyte_mass : float
        Mass of analyte in kg.
    titrant_molinity : float
        Molinity of titrant in mol/kg-sol.
    titrant_normality : float, optional
        Titrant normality, by default 1 (e.g., for HCl).

    Returns
    -------
    GranGuessesResult : namedtuple with the fields
        alkalinity : float
            Alkalinity estimate in mol/kg-sol.
        gfunc : array-like float
            Gran function values from `gran_function`.
        used : array-like bool
            Which Gran function points are used.
        intercept_x : float
            x-axis intercept computed from `lr`.
        lr : scipy.stats.LinregressResult
            Output from `scipy.stats.linregress`.
        emf0 : float
            Final Gran-plot estimate of EMF0 in mV.
        emf0s : float
            Gran-plot estimates of EMF0 in mV following DAA03 eq. 11.
        pH : array-like (float)
            pH through the titration based on estimated EMF0.
    """
    ga = gran_alkalinity(
        titrant_mass,
        emf,
        temperature,
        analyte_mass,
        titrant_molinity,
        titrant_normality=titrant_normality,
    )
    emf0s = gran_emf0s(
        titrant_mass,
        emf,
        temperature,
        analyte_mass,
        titrant_molinity,
        ga.alkalinity,
        titrant_normality=titrant_normality,
    )
    emf0 = np.mean(emf0s[ga.used])
    pH = convert.emf_to_pH(emf, emf0, temperature)
    return GranGuessesResult(*ga, emf0, emf0s, pH)


def totals_ks(converted, **kwargs):
    """Get total salt contents and equilibrium constants through a titration.

    Parameters
    ----------
    converted : Converted
        A namedtuple generated by `convert.amount_units`.
    salinity : float
        Practical salinity of the analyte at the start of the titration.
    kwargs - additional kwargs for PyCO2SYS (see its docs for detail):
        dic
        k_alpha
        k_ammonia
        k_beta
        k_bisulfate
        k_borate
        k_carbonic_1
        k_carbonic_2
        k_fluoride
        k_phosphoric_1
        k_phosphoric_2
        k_phosphoric_3
        k_silicate
        k_sulfide
        k_water
        opt_k_bisulfate
        opt_k_carbonic
        opt_k_fluoride
        opt_pH_scale
        opt_total_borate
        total_alpha
        total_ammonia
        total_beta
        total_borate
        total_fluoride
        total_phosphate
        total_silicate
        total_sulfate
        total_sulfide

    Returns
    -------
    totals : dict
        The total salt contents through the titration, including dilution
        by the titrant, assuming the titrant is not one of the totals (e.g.,
        it is not H2SO4).
    k_constants : dict
        The equilibrium constants through the titration, including dilution
        by the titrant (affects pH scale conversions only) and temperature
        variations.
    """
    cv = converted
    # Get totals from PyCO2SYS
    kwargs_totals = {
        k: v
        for k, v in kwargs.items()
        if k
        in [
            "dic",
            "total_alpha",
            "total_beta",
            "total_ammonia",
            "total_phosphate",
            "total_silicate",
            "total_sulfide",
            "total_borate",
            "total_fluoride",
            "total_sulfate",
            "opt_k_carbonic",
            "opt_total_borate",
        ]
    }
    totals, totals_pyco2 = interface.get_totals(cv.salinity, **kwargs_totals)
    # Dilute totals with titrant
    totals = convert.dilute_totals(totals, cv.titrant_mass, cv.analyte_mass)
    totals_pyco2 = convert.dilute_totals_pyco2(
        totals_pyco2, cv.titrant_mass, cv.analyte_mass
    )
    # Get k_constants from PyCO2SYS
    kwargs_k_constants = {
        k: v
        for k, v in kwargs.items()
        if k
        in [
            "k_alpha",
            "k_ammonia",
            "k_beta",
            "k_bisulfate",
            "k_borate",
            "k_carbonic_1",
            "k_carbonic_2",
            "k_fluoride",
            "k_phosphoric_1",
            "k_phosphoric_2",
            "k_phosphoric_3",
            "k_silicate",
            "k_sulfide",
            "k_water",
            "opt_k_bisulfate",
            "opt_k_carbonic",
            "opt_k_fluoride",
            "opt_pH_scale",
            "opt_total_borate",
        ]
    }
    k_constants = interface.get_k_constants(
        totals_pyco2, cv.temperature, **kwargs_k_constants
    )
    return totals, k_constants


def add_titrant_totals(
    totals,
    titrant_mass,
    analyte_mass,
    titrant_molinity,
    titrant_molinity_prev=0,
    inplace=False,
    **titrant_totals,
):
    if not inplace:
        totals = totals.copy()
    for t_total, factor in titrant_totals.items():
        assert t_total.startswith("titrant_total_")
        total = t_total[8:]
        totals[total] += (
            titrant_mass
            * (titrant_molinity - titrant_molinity_prev)
            * factor
            / (titrant_mass + analyte_mass)
        )
    if not inplace:
        return totals


def _lsqfun_solve_emf(
    alkalinity_emf0,
    titrant_molinity,
    titrant_mass,
    emf,
    temperature,
    analyte_mass,
    totals,
    k_constants,
    titrant_normality,
):
    """Calculate residuals for the solver."""
    alkalinity, emf0 = alkalinity_emf0
    pH = convert.emf_to_pH(emf, emf0, temperature)
    mixture_mass = titrant_mass + analyte_mass
    dilution_factor = convert.get_dilution_factor(titrant_mass, analyte_mass)
    return (
        simulate.alkalinity(pH, totals, k_constants)
        - alkalinity * dilution_factor
        + titrant_mass * titrant_molinity * titrant_normality / mixture_mass
    )


def solve_emf(
    titrant_molinity,
    titrant_mass,
    emf,
    temperature,
    analyte_mass,
    totals,
    k_constants,
    alkalinity_initial=None,
    emf0_initial=None,
    pH_min=3,
    pH_max=4,
    titrant_normality=1,
):
    """Solve for alkalinity and EMF0 using the complete-calculation method,
    assuming a titrant normality of 1 (e.g., HCl), when EMF is known.

    Parameters
    ----------
    titrant_molinity : float
        Molinity of titrant in mol/kg-sol.
    titrant_mass : array-like float
        Mass of titrant in kg.
    emf : array-like float
        EMF measured across the titrant-analyte mixture in mV.
    temperature : array-like float
        Temperature of titrant-analyte mixture in °C.
    analyte_mass : float
        Mass of analyte in kg.
    totals : dict of array-like floats
        Total salt contents through the titration, created with `totals_ks` and
        with any additions from the titrant included with `add_titrant_totals`.
    k_constants : dict of array-like floats
        Equilibrium constants through the titration, created with `totals_ks`.
    alkalinity_initial : float, optional
        An alkalinity value in µmol/kg-sol to use to initialise the
        solver.  By default None, in which case this is estimated using the
        Gran approach (see `gran_guesses`).
    emf0_initial : float, optional
        An EMF0 value to use to calculate the initial pH estimates from EMF,
        which are used with pH_min and pH_max to find the data points to use
        for solving.  By default None, in which case this is estimated using
        the Gran approach (see `gran_guesses`).
    pH_min : float, optional
        Minimum pH to use from the titration data, by default 3.
    pH_max : float, optional
        Maximum pH to use from the titration data, by default 4.
    titrant_normality : float, optional
        Titrant normality, by default 1 (e.g., for HCl).

    Returns
    -------
    SolveEmfResult : namedtuple with the fields
        alkalinity : float
            Total alkalinity in µmol/kg-sol.
        emf0 : float
            EMF0 in mV.
        used : array-like bool
            Which data points were used.
        pH : array-like float
            pH on the free scale.
        alkalinity_std : float
            Standard deviation of alkalinity estimates in µmol/kg-sol.
        alkalinity_all : array-like float
            Alkalinity estimates at every titration point in µmol/kg-sol.
        opt_result : scipy.optimize.OptimizeResult
            Output from `scipy.optimize.least_squares`.
        ggr : GranGuessesResult
            Output from `gran_guesses`.
    """
    # Get initial guesses
    ggr = gran_guesses(
        titrant_mass,
        emf,
        temperature,
        analyte_mass,
        titrant_molinity,
        titrant_normality=titrant_normality,
    )
    if alkalinity_initial is None:
        alkalinity = ggr.alkalinity
    else:
        alkalinity = alkalinity_initial * 1e-6
    # Set which data points to use in the final solver
    assert pH_min < pH_max
    if emf0_initial is None:
        emf0 = ggr.emf0
        pH = ggr.pH
    else:
        emf0 = emf0_initial
        pH = convert.emf_to_pH(emf, emf0, temperature)
    used = (pH >= pH_min) & (pH <= pH_max)
    totals_used = {
        k: v[used] if np.size(v) > 1 else v for k, v in totals.items()
    }
    ks_used = {
        k: v[used] if np.size(v) > 1 else v for k, v in k_constants.items()
    }
    # Solve for alkalinity and EMF0
    opt_result = least_squares(
        _lsqfun_solve_emf,
        [alkalinity, emf0],
        args=(
            titrant_molinity,
            titrant_mass[used],
            emf[used],
            temperature[used],
            analyte_mass,
            totals_used,
            ks_used,
            titrant_normality,
        ),
        **kwargs_least_squares,
    )
    # Unpack and process results
    alkalinity = opt_result["x"][0] * 1e6
    alkalinity_std = np.std(opt_result["fun"]) * 1e6
    emf0 = opt_result["x"][1]
    pH = convert.emf_to_pH(emf, emf0, temperature)
    alkalinity_all = (
        1e6
        * (
            simulate.alkalinity(pH, totals, k_constants)
            + (titrant_mass * titrant_molinity * titrant_normality)
            / (titrant_mass + analyte_mass)
        )
        / convert.get_dilution_factor(titrant_mass, analyte_mass)
    )
    return SolveEmfResult(
        alkalinity,
        emf0,
        used,
        pH,
        alkalinity_std,
        alkalinity_all,
        opt_result,
        ggr,
    )


def solve_pH(
    titrant_molinity,
    titrant_mass,
    pH,
    temperature,
    analyte_mass,
    totals,
    k_constants,
    pH_min=3,
    pH_max=4,
    titrant_normality=1,
):
    """Solve for alkalinity and EMF0 using the complete-calculation method,
    assuming a titrant normality of 1 (e.g., HCl), when pH is known.

    Parameters
    ----------
    titrant_molinity : float
        Molinity of titrant in mol/kg-sol.
    titrant_mass : array-like float
        Mass of titrant in kg.
    pH : array-like float
        pH in the titrant-analyte mixture on the same scale as `k_constants`.
    temperature : array-like float
        Temperature of titrant-analyte mixture in °C.
    analyte_mass : float
        Mass of analyte in kg.
    totals : dict of array-like floats
        Total salt contents through the titration, created with `totals_ks` and
        with any additions from the titrant included with `add_titrant_totals`.
    k_constants : dict of array-like floats
        Equilibrium constants through the titration, created with `totals_ks`,
        and on the same scale as the `pH`.
    pH_min : float, optional
        Minimum pH to use from the titration data, by default 3.
    pH_max : float, optional
        Maximum pH to use from the titration data, by default 4.
    titrant_normality : float, optional
        Titrant normality, by default 1 (e.g., for HCl).

    Returns
    -------
    SolvePhResult : namedtuple with the fields
        alkalinity : float
            Total alkalinity in µmol/kg-sol.
        used : array-like bool
            Which data points were used.
        alkalinity_std : float
            Standard deviation of alkalinity estimates in µmol/kg-sol.
        alkalinity_all : array-like float
            Alkalinity estimates at every titration point in µmol/kg-sol.
    """
    assert pH_min < pH_max
    used = (pH >= pH_min) & (pH <= pH_max)
    alkalinity_all = (
        1e6
        * (
            simulate.alkalinity(pH, totals, k_constants)
            + (titrant_mass * titrant_molinity * titrant_normality)
            / (titrant_mass + analyte_mass)
        )
        / convert.get_dilution_factor(titrant_mass, analyte_mass)
    )
    alkalinity = np.mean(alkalinity_all[used])
    alkalinity_std = np.std(alkalinity_all[used])
    return SolvePhResult(alkalinity, used, alkalinity_std, alkalinity_all)


def _lsqfun_calibrate_emf(
    titrant_molinity,
    alkalinity_certified,
    titrant_mass,
    measurement,
    temperature,
    analyte_mass,
    totals,
    k_constants,
    alkalinity_initial,
    emf0_initial,
    pH_min,
    pH_max,
    titrant_normality,
    **titrant_totals,
):
    """Calculate residuals for the calibrator."""
    # Add titrant to totals (only relevant for H2SO4 etc. titrant)
    totals = add_titrant_totals(
        totals,
        titrant_mass,
        analyte_mass,
        titrant_molinity,
        titrant_molinity_prev=0,
        **titrant_totals,
    )
    # Solve for alkalinity and EMF
    sr = solve_emf(
        titrant_molinity[0],
        titrant_mass,
        measurement,
        temperature,
        analyte_mass,
        totals,
        k_constants,
        alkalinity_initial=alkalinity_initial,
        emf0_initial=emf0_initial,
        pH_min=pH_min,
        pH_max=pH_max,
        titrant_normality=titrant_normality,
    )
    # Revert to original totals
    totals = add_titrant_totals(
        totals,
        titrant_mass,
        analyte_mass,
        0,
        titrant_molinity_prev=titrant_molinity,
        **titrant_totals,
    )
    return sr.alkalinity - alkalinity_certified


def calibrate_emf(
    alkalinity_certified,
    titrant_mass,
    emf,
    temperature,
    analyte_mass,
    totals,
    k_constants,
    alkalinity_initial=None,
    emf0_initial=None,
    pH_min=3,
    pH_max=4,
    titrant_molinity_initial=0.1,
    titrant_normality=1,
    **titrant_totals,
):
    """Solve for `titrant_molinity` given `alkalinity_certified`.

    Parameters
    ----------
    alkalinity_certified : float
        The target total alkalinity in µmol/kg-sol.
    titrant_mass : array-like
        Mass of titrant in kg.
    emf : array-like
        EMF measured across the titrant-analyte mixture in mV.
    temperature : array-like
        Temperature of titrant-analyte mixture in °C.
    analyte_mass : float
        Mass of analyte in kg.
    totals : dict
        Total salt contents through the titration, created with
        `totals_ks`.  Any additional totals added by the titrant should
        NOT have been added here yet (unlike for `solve`).
    k_constants : dict
        Equilibrium constants through the titration, created with
        `totals_ks`.
    alkalinity_initial : float, optional
        An alkalinity value in µmol/kg-sol to use to initialise the
        solver.  By default None, in which case this is estimated using the
        Gran approach (see `gran_guesses`).
    emf0_initial : float, optional
        An EMF0 value to use to calculate the initial pH estimates from EMF,
        which are used with pH_min and pH_max to find the data points to use
        for solving.  By default None, in which case this is estimated using
        the Gran approach (see `gran_guesses`).
    pH_min : float, optional
        Minimum pH to use from the titration data, by default 3.
    pH_max : float, optional
        Maximum pH to use from the titration data, by default 4.
    titrant_molinity_initial : float, optional
        First guess for the molinity of titrant in mol/kg-sol, by default
        0.1.
    titrant_normality : float, optional
        Titrant normality, by default 1 (e.g., for HCl).

    Returns
    -------
    opt_result : scipy.optimize.OptimizeResult
        Output from `scipy.optimize.least_squares`, where the solved value is
        `titrant_molinity = opt_result["x"][0]`.
    """
    return least_squares(
        _lsqfun_calibrate_emf,
        [titrant_molinity_initial],
        args=(
            alkalinity_certified,
            titrant_mass,
            emf,
            temperature,
            analyte_mass,
            totals,
            k_constants,
            alkalinity_initial,
            emf0_initial,
            pH_min,
            pH_max,
            titrant_normality,
        ),
        kwargs=titrant_totals,
        **kwargs_least_squares,
    )


def _lsqfun_calibrate_pH(
    titrant_molinity,
    alkalinity_certified,
    titrant_mass,
    pH,
    temperature,
    analyte_mass,
    totals,
    k_constants,
    pH_min,
    pH_max,
    titrant_normality,
    **titrant_totals,
):
    """Calculate residuals for the calibrator."""
    # Add titrant to totals (only relevant for H2SO4 etc. titrant)
    totals = add_titrant_totals(
        totals,
        titrant_mass,
        analyte_mass,
        titrant_molinity,
        titrant_molinity_prev=0,
        **titrant_totals,
    )
    # Solve for alkalinity
    sr = solve_pH(
        titrant_molinity[0],
        titrant_mass,
        pH,
        temperature,
        analyte_mass,
        totals,
        k_constants,
        pH_min=pH_min,
        pH_max=pH_max,
        titrant_normality=titrant_normality,
    )
    # Revert to original totals
    totals = add_titrant_totals(
        totals,
        titrant_mass,
        analyte_mass,
        0,
        titrant_molinity_prev=titrant_molinity,
        **titrant_totals,
    )
    return sr.alkalinity - alkalinity_certified


def calibrate_pH(
    alkalinity_certified,
    titrant_mass,
    pH,
    temperature,
    analyte_mass,
    totals,
    k_constants,
    pH_min=3,
    pH_max=4,
    titrant_molinity_initial=0.1,
    titrant_normality=1,
    **titrant_totals,
):
    """Solve for `titrant_molinity` given `alkalinity_certified`.

    Parameters
    ----------
    alkalinity_certified : float
        The target total alkalinity in µmol/kg-sol.
    titrant_mass : array-like
        Mass of titrant in kg.
    pH : array-like
        pH in the titrant-analyte mixture on the same scale as `k_constants`.
    temperature : array-like
        Temperature of titrant-analyte mixture in °C.
    analyte_mass : float
        Mass of analyte in kg.
    totals : dict
        Total salt contents through the titration, created with
        `totals_ks`.  Any additional totals added by the titrant should
        NOT have been added here yet (unlike for `solve`).
    k_constants : dict
        Equilibrium constants through the titration, created with
        `totals_ks`, and on the same pH scale as the `pH`.
    pH_min : float, optional
        Minimum pH to use from the titration data, by default 3.
    pH_max : float, optional
        Maximum pH to use from the titration data, by default 4.
    titrant_molinity_initial : float, optional
        First guess for the molinity of titrant in mol/kg-sol, by default
        0.1.
    titrant_normality : float, optional
        Titrant normality, by default 1 (e.g., for HCl).

    Returns
    -------
    opt_result : scipy.optimize.OptimizeResult
        Output from `scipy.optimize.least_squares`, where the solved value is
        `titrant_molinity = opt_result["x"][0]`.
    """
    return least_squares(
        _lsqfun_calibrate_pH,
        [titrant_molinity_initial],
        args=(
            alkalinity_certified,
            titrant_mass,
            pH,
            temperature,
            analyte_mass,
            totals,
            k_constants,
            pH_min,
            pH_max,
            titrant_normality,
        ),
        kwargs=titrant_totals,
        **kwargs_least_squares,
    )
