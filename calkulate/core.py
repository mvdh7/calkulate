# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019--2025  Matthew P. Humphreys  (GNU GPLv3)
"""Calibrate and solve titration datasets.

This is the lowest level of the three layers of processing functions.
Titration data need to be imported and converted into separate NumPy arrays in
order to work with the functions here.

First-guess functions
---------------------
gran_estimator
gran_guess_alkalinity
gran_guesses_emf0
gran_guesses

Solver functions
----------------
solve_emf_complete
solve_emf_complete_H2SO4
solve_emf_pH_adjust

Calibration functions
---------------------
calibrate
"""

from collections import namedtuple

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import linregress

from . import constants, convert, simulate
from .settings import least_squares_kwargs


GranGuessesResult = namedtuple(
    "GranGuessesResult", ("alkalinity", "emf0", "pH", "G")
)


def gran_estimator(mixture_mass, emf, temperature):
    """Calculate Gran-plot estimator (DAA03 eq. 10) using all provided data.

    Parameters
    ----------
    mixture_mass : float
        Mass of titrant-analyte mixture in kg.
    emf : float
        EMF measured in titrant-analyte mixture in mV.
    temperature : float
        Temperature of titrant-analyte mixture in °C.

    Returns
    -------
    float
        Gran-plot estimator (DAA03 eq. 10).
    """
    return mixture_mass * np.exp(
        emf
        * constants.faraday
        / (constants.ideal_gas * (temperature + constants.absolute_zero))
    )


def gran_guess_alkalinity(
    titrant_mass,
    gran_estimates,
    analyte_mass,
    titrant_molinity,
    titrant_normality=1,
):
    """Gran-plot first guess of alkalinity using all provided data.

    Parameters
    ----------
    titrant_mass : array-like
        Mass of titrant in kg.
    gran_estimates : array-like
        Gran-plot estimator, output from `calkulate.core.gran_estimator`.
    analyte_mass : float
        Mass of analyte in kg.
    titrant_molinity : float
        Molinity of titrant in mol/kg-solution.
    titrant_normality : float, optional
        Titrant normality, by default 1 (e.g., for HCl).

    Returns
    -------
    alkalinity_guess : float
        Gran-plot estimate of alkalinity in mol/kg-solution.
    """

    # Do regression through simple Gran-plot estimator
    lr = linregress(titrant_mass, gran_estimates)
    # Find alkalinity guess from the x-axis intercept
    intercept_x = -lr.intercept / lr.slope
    alkalinity_guess = (
        intercept_x * titrant_molinity * titrant_normality / analyte_mass
    )
    return alkalinity_guess


def gran_guesses_emf0(
    titrant_mass,
    emf,
    temperature,
    analyte_mass,
    titrant_molinity,
    alkalinity_guess=None,
    titrant_normality=1,
    total_fluoride=0,
    total_sulfate=0,
):
    """Gran-plot first guesses of EMF0 (DAA03 eq. 11) using all provided data.

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
        Molinity of titrant in mol/kg-solution.
    alkalinity_guess : float, optional
        First-guess value for total alkalinity in mol/kg-solution, by default
        `None`, in which case `gran_guess_alkalinity` is used.
    titrant_normality : float, optional
        Titrant normality, by default 1 (e.g., for HCl).
    total_fluoride : float, optional
        Total fluoride in titrant-analyte mixture in mol/kg-solution, by
        default 0.
    total_sulfate : float, optional
        Total sulfate in titrant-analyte mixture in mol/kg-solution, by
        default 0.

    Returns
    -------
    emf0_guesses : array-like or float
        Gran-plot guesses of EMF0 in mV following DAA03 eq. 11.
    """
    # Get alkalinity_guess if one is not already provided
    mixture_mass = titrant_mass + analyte_mass
    if alkalinity_guess is None:
        gran_estimates = gran_estimator(mixture_mass, emf, temperature)
        alkalinity_guess = gran_guess_alkalinity(
            titrant_mass,
            gran_estimates,
            analyte_mass,
            titrant_molinity,
            titrant_normality=titrant_normality,
        )
    # Calculate first-guesses of EMF0
    temperature_K = temperature + constants.absolute_zero
    emf0_guesses = emf - (
        constants.ideal_gas * temperature_K / constants.faraday
    ) * np.log(
        (
            (
                titrant_mass * titrant_molinity * titrant_normality
                - analyte_mass * alkalinity_guess
            )
            - analyte_mass * (total_fluoride + total_sulfate)
        )
        / mixture_mass
    )
    return emf0_guesses


def gran_guesses(
    titrant_mass,
    emf,
    temperature,
    analyte_mass,
    titrant_molinity,
    emf0_guess=None,
    titrant_normality=1,
    total_fluoride=0,
    total_sulfate=0,
):
    """Calculate Gran-plot first guesses for alkalinity, EMF0 and pH, using a
    subset of the provided data points.

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
        Molinity of titrant in mol/kg-solution.
    emf0_guess : float, optional
        Gran-plot guess of EMF0 in mV, by default `None`, in which case it is
        determined using `gran_guesses_emf0`.
    titrant_normality : float, optional
        Titrant normality, by default 1 (e.g., for HCl).
    total_fluoride : float, optional
        Total fluoride in titrant-analyte mixture in mol/kg-solution, by
        default 0.
    total_sulfate : float, optional
        Total sulfate in titrant-analyte mixture in mol/kg-solution, by
        default 0.

    Returns
    -------
    GranGuessesResult
        A namedtuple containing the fields
            alkalinity : float
                Gran-plot estimate of alkalinity in mol/kg-solution.
            emf0 : float
                Gran-plot estimate of EMF0 in mV following DAA03 eq. 11.
            pH : array-like (float)
                pH through the titration based on estimated EMF0.
            G : array-like (bool)
                Which titration data points were used to estimate alkalinity
                and EMF0.
    """
    # Get simple Gran-plot estimator
    mixture_mass = titrant_mass + analyte_mass
    gran_estimates = gran_estimator(mixture_mass, emf, temperature)
    # Select which data points to use for first guesses
    G = (gran_estimates >= 0.1 * np.max(gran_estimates)) & (
        gran_estimates <= 0.9 * np.max(gran_estimates)
    )
    # Make first guesses
    alkalinity_guess = gran_guess_alkalinity(
        titrant_mass[G],
        gran_estimates[G],
        analyte_mass,
        titrant_molinity,
        titrant_normality=titrant_normality,
    )
    if np.size(temperature) == 1:
        temperature = np.full(np.size(titrant_mass), temperature)
    if emf0_guess is None:
        emf0_guess = np.mean(
            gran_guesses_emf0(
                titrant_mass[G],
                emf[G],
                temperature[G],
                analyte_mass,
                titrant_molinity,
                alkalinity_guess=alkalinity_guess,
                titrant_normality=titrant_normality,
                total_fluoride=total_fluoride,
                total_sulfate=total_sulfate,
            )
        )
    pH_guesses = convert.emf_to_pH(emf, emf0_guess, temperature)
    return GranGuessesResult(alkalinity_guess, emf0_guess, pH_guesses, G)


def _lsqfun_solve_emf_complete(
    alkalinity_emf0,
    titrant_molinity,
    titrant_mass,
    emf,
    temperature,
    analyte_mass,
    totals,
    k_constants,
):
    """Calculate residuals for the solvers `solve_emf_complete` and
    `solve_emf_pH_adjust`.
    """
    alkalinity, emf0 = alkalinity_emf0
    pH = convert.emf_to_pH(emf, emf0, temperature)
    mixture_mass = titrant_mass + analyte_mass
    dilution_factor = convert.get_dilution_factor(titrant_mass, analyte_mass)
    return (
        simulate.alkalinity(pH, totals, k_constants)
        - alkalinity * dilution_factor
        + titrant_mass * titrant_molinity / mixture_mass
    )


def solve_emf_complete(
    titrant_molinity,
    titrant_mass,
    emf,
    temperature,
    analyte_mass,
    totals,
    k_constants,
    emf0_guess=None,
    pH_min=3,
    pH_max=4,
):
    """Solve for alkalinity and EMF0 using the complete-calculation method,
    assuming a titrant normality of 1 (e.g., HCl).

    Parameters
    ----------
    titrant_molinity : float
        Molinity of titrant in mol/kg-solution.
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
        `titration.get_totals_k_constants`.
    k_constants : dict
        Equilibrium constants through the titration, created with
        `titration.get_totals_k_constants`.
    emf0_guess : float, optional
        Gran-plot guess of EMF0 in mV, by default `None`, in which case it is
        determined using `gran_guesses_emf0`.
    pH_min : float, optional
        Minimum pH to use from the titration data, by default 3.
    pH_max : float, optional
        Maximum pH to use from the titration data, by default 4.

    Returns
    -------
    opt_result : dict
        The output from `scipy.optimize.least_squares`, where solved values
        are `alkalinity, emf0 = opt_result["x"]`, plus plus some additional
        entries:
            data_used : bool
                Logical for which titration data points were used (i.e., in
                the range from `pH_min` to `pH_max`).
            alkalinity_guess : float
                First-guess alkalinity value in µmol/kg-solution.
            emf0_guess : float
                First-guess EMF0 value in µmol/kg-solution.
            emf0_offset : float
                Not used here (contains `np.nan`).
    """
    # Get initial guesses
    ggr = gran_guesses(
        titrant_mass,
        emf,
        temperature,
        analyte_mass,
        titrant_molinity,
        emf0_guess=emf0_guess,
    )
    # Set which data points to use in the final solver
    assert pH_min < pH_max
    G = (ggr.pH >= pH_min) & (ggr.pH <= pH_max)
    totals_G = {k: v[G] if np.size(v) > 1 else v for k, v in totals.items()}
    k_constants_G = {
        k: v[G] if np.size(v) > 1 else v for k, v in k_constants.items()
    }
    # Solve for alkalinity and EMF0
    opt_result = least_squares(
        _lsqfun_solve_emf_complete,
        [ggr.alkalinity, ggr.emf0],
        args=(
            titrant_molinity,
            titrant_mass[G],
            emf[G],
            temperature[G],
            analyte_mass,
            totals_G,
            k_constants_G,
        ),
        x_scale=[1e-6, 1],
        **least_squares_kwargs,
    )
    # Append which data points were used and initial guesses to the output
    opt_result["data_used"] = G
    opt_result["alkalinity_guess"] = ggr.alkalinity * 1e6
    opt_result["emf0_guess"] = ggr.emf0
    opt_result["emf0_offset"] = np.nan
    return opt_result


def solve_emf_pH_adjust(
    titrant_molinity,
    titrant_mass,
    pH,
    temperature,
    analyte_mass,
    totals,
    k_constants,
    emf0_guess=None,
    pH_min=3,
    pH_max=4,
):
    """Solve for alkalinity and ∆EMF0 when pH is known but may be adjusted,
    assuming a titrant normality of 1 (e.g., HCl).

    Parameters
    ----------
    titrant_molinity : float
        Molinity of titrant in mol/kg-solution.
    titrant_mass : array-like
        Mass of titrant in kg.
    pH : array-like
        pH in the titrant-analyte mixture through the titration.
    temperature : array-like
        Temperature of titrant-analyte mixture in °C.
    analyte_mass : float
        Mass of analyte in kg.
    totals : dict
        Total salt contents through the titration, created with
        `titration.get_totals_k_constants`.
    k_constants : dict
        Equilibrium constants through the titration, created with
        `titration.get_totals_k_constants`.
    emf0_guess : float, optional
        Not used here.
    pH_min : float, optional
        Minimum pH to use from the titration data, by default 3.
    pH_max : float, optional
        Maximum pH to use from the titration data, by default 4.

    Returns
    -------
    opt_result : dict
        The output from `scipy.optimize.least_squares`, where solved values
        are `alkalinity, emf0_offset = opt_result["x"]`, plus plus some
        additional entries:
            data_used : bool
                Logical for which titration data points were used (i.e., in
                the range from `pH_min` to `pH_max`).
            alkalinity_guess : float
                First-guess alkalinity value in µmol/kg-solution.
            emf0_guess : float
                Not used here (contains `np.nan`).
            emf0_offset : float
                Adjustment applied to EMF0 in mV.
    """
    # Set which data points to use in the final solver
    assert pH_min < pH_max
    G = (pH >= pH_min) & (pH <= pH_max)
    totals_G = {k: v[G] if np.size(v) > 1 else v for k, v in totals.items()}
    k_constants_G = {
        k: v[G] if np.size(v) > 1 else v for k, v in k_constants.items()
    }
    alkalinity_guess = (
        (
            simulate.alkalinity(pH[G], totals_G, k_constants_G)
            + titrant_mass[G]
            * titrant_molinity
            / (titrant_mass[G] + analyte_mass)
        )
        / convert.get_dilution_factor(titrant_mass[G], analyte_mass)
    ).mean()
    # Solve for alkalinity and EMF0
    opt_result = least_squares(
        _lsqfun_solve_emf_complete,
        [alkalinity_guess, 0],
        args=(
            titrant_molinity,
            titrant_mass[G],
            convert.pH_to_emf(pH, 0, temperature)[G],
            temperature[G],
            analyte_mass,
            totals_G,
            k_constants_G,
        ),
        x_scale=[1e-6, 1],
        **least_squares_kwargs,
    )
    # Add which data points were used and initial guesses to the output
    opt_result["data_used"] = G
    opt_result["alkalinity_guess"] = alkalinity_guess * 1e6
    opt_result["emf0_guess"] = np.nan
    opt_result["emf0_offset"] = 0
    return opt_result


def _lsqfun_solve_emf_complete_H2SO4(
    alkalinity_emf0,
    titrant_molinity,
    titrant_mass,
    emf,
    temperature,
    analyte_mass,
    totals,
    k_constants,
):
    """Calculate residuals for `solve_emf_complete_H2SO4`."""
    alkalinity, emf0 = alkalinity_emf0
    pH = convert.emf_to_pH(emf, emf0, temperature)
    mixture_mass = titrant_mass + analyte_mass
    dilution_factor = convert.get_dilution_factor(titrant_mass, analyte_mass)
    residual = (
        simulate.alkalinity(pH, totals, k_constants)
        - alkalinity * dilution_factor
        + 2 * titrant_mass * titrant_molinity / mixture_mass
    )
    return residual


def solve_emf_complete_H2SO4(
    titrant_molinity,
    titrant_mass,
    emf,
    temperature,
    analyte_mass,
    totals,
    k_constants,
    emf0_guess=None,
    pH_min=3,
    pH_max=4,
):
    """Solve for alkalinity and EMF0 using the complete-calculation method
    for an H2SO4 titrant.


    Parameters
    ----------
    titrant_molinity : float
        Molinity of titrant in mol/kg-solution.
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
        `titration.get_totals_k_constants`.
    k_constants : dict
        Equilibrium constants through the titration, created with
        `titration.get_totals_k_constants`.
    emf0_guess : float, optional
        Gran-plot guess of EMF0 in mV, by default `None`, in which case it is
        determined using `gran_guesses_emf0`.
    pH_min : float, optional
        Minimum pH to use from the titration data, by default 3.
    pH_max : float, optional
        Maximum pH to use from the titration data, by default 4.

    Returns
    -------
    opt_result : dict
        The output from `scipy.optimize.least_squares`, where solved values
        are `alkalinity, emf0 = opt_result["x"]`, plus plus some additional
        entries:
            data_used : bool
                Logical for which titration data points were used (i.e., in
                the range from `pH_min` to `pH_max`).
            alkalinity_guess : float
                First-guess alkalinity value in µmol/kg-solution.
            emf0_guess : float
                First-guess EMF0 value in µmol/kg-solution.
            emf0_offset : float
                Not used here (contains `np.nan`).
    """
    # Get initial guesses
    ggr = gran_guesses(
        titrant_mass,
        emf,
        temperature,
        analyte_mass,
        titrant_molinity,
        emf0_guess=emf0_guess,
        titrant_normality=2,
    )[:-1]
    # Set which data points to use in the final solver
    assert pH_min < pH_max
    G = (ggr.pH >= pH_min) & (ggr.pH <= pH_max)
    totals_G = {k: v[G] if np.size(v) > 1 else v for k, v in totals.items()}
    k_constants_G = {
        k: v[G] if np.size(v) > 1 else v for k, v in k_constants.items()
    }
    # Solve for alkalinity and EMF0
    opt_result = least_squares(
        _lsqfun_solve_emf_complete_H2SO4,
        [ggr.alkalinity, ggr.emf0],
        args=(
            titrant_molinity,
            titrant_mass[G],
            emf[G],
            temperature[G],
            analyte_mass,
            totals_G,
            k_constants_G,
        ),
        x_scale=[1e-6, 1],
        **least_squares_kwargs,
    )
    # Add which data points were used to the output
    opt_result["data_used"] = G
    opt_result["alkalinity_guess"] = ggr.alkalinity * 1e6
    opt_result["emf0_guess"] = ggr.emf0
    opt_result["emf0_offset"] = np.nan
    return opt_result


def _lsqfun_calibrate(
    titrant_molinity,
    alkalinity_certified,
    titrant_mass,
    measurement,
    temperature,
    analyte_mass,
    totals,
    k_constants,
    pH_min,
    pH_max,
    solve_func,
):
    """Calculate residuals for the calibrator."""
    alkalinity = solve_func(
        titrant_molinity[0],
        titrant_mass,
        measurement,
        temperature,
        analyte_mass,
        totals,
        k_constants,
        pH_min=pH_min,
        pH_max=pH_max,
    )["x"][0]
    return alkalinity * 1e6 - alkalinity_certified


def calibrate(
    alkalinity_certified,
    titrant_mass,
    measurement,
    temperature,
    analyte_mass,
    totals,
    k_constants,
    pH_min=3,
    pH_max=4,
    solve_func=solve_emf_complete,
    titrant_molinity_guess=0.1,
):
    """Solve for `titrant_molinity` given `alkalinity_certified`.


    Parameters
    ----------
    alkalinity_certified : float
        The target total alkalinity in µmol/kg-solution.
    titrant_mass : array-like
        Mass of titrant in kg.
    measurement : array-like
        Depending on what is required by `solve_func`, either EMF measured
        across the titrant-analyte mixture in mV, or its pH.
    temperature : array-like
        Temperature of titrant-analyte mixture in °C.
    analyte_mass : float
        Mass of analyte in kg.
    totals : dict
        Total salt contents through the titration, created with
        `titration.get_totals_k_constants`.
    k_constants : dict
        Equilibrium constants through the titration, created with
        `titration.get_totals_k_constants`.
    pH_min : float, optional
        Minimum pH to use from the titration data, by default 3.
    pH_max : float, optional
        Maximum pH to use from the titration data, by default 4.
    solve_func : func, optional
        Which function to use to solve for alkalinity, by default
        `solve_emf_complete`.
    titrant_molinity_guess : float, optional
        First guess for the molinity of titrant in mol/kg-solution, by default
        0.1.

    Returns
    -------
    opt_result : dict
        The output from `scipy.optimize.least_squares`, where solved values
        are `titrant_molinity = opt_result["x"][0]`.
    """
    return least_squares(
        _lsqfun_calibrate,
        [titrant_molinity_guess],
        args=(
            alkalinity_certified,
            titrant_mass,
            measurement,
            temperature,
            analyte_mass,
            totals,
            k_constants,
            pH_min,
            pH_max,
            solve_func,
        ),
        **least_squares_kwargs,
    )
