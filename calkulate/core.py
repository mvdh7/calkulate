# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019--2021  Matthew P. Humphreys  (GNU GPLv3)
"""Calibrate and solve titration datasets."""

import copy
import numpy as np
from scipy.stats import linregress
from scipy.optimize import least_squares
from . import constants, convert, default, interface, simulate


def gran_estimator(mixture_mass, emf, temperature):
    """Simple Gran-plot estimator following DAA03 eq. 10.
    Uses all provided data points.
    """
    gran_estimates = mixture_mass * np.exp(
        emf
        * constants.faraday
        / (constants.ideal_gas * (temperature + constants.absolute_zero))
    )
    return gran_estimates


def gran_guess_alkalinity(
    titrant_mass,
    gran_estimates,
    analyte_mass,
    titrant_molinity,
    titrant=default.titrant,
):
    """Simple Gran-plot first-guess of alkalinity.
    Uses all provided data points.
    """
    # Do regression through simple Gran-plot estimator
    gradient, intercept_y = linregress(titrant_mass, gran_estimates)[:2]
    # Find alkalinity guess from the x-axis intercept
    intercept_x = -intercept_y / gradient
    alkalinity_guess = intercept_x * titrant_molinity / analyte_mass
    return alkalinity_guess, gradient, intercept_y


def gran_guesses_emf0(
    titrant_mass,
    emf,
    temperature,
    analyte_mass,
    titrant_molinity,
    alkalinity_guess=None,
    titrant=default.titrant,
    HF=0,
    HSO4=0,
):
    """Simple Gran-plot first-guess of EMF0 following DAA03 eq. 11.
    Uses all provided data points.
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
        )[0]
        if titrant == "H2SO4":
            alkalinity_guess *= 2
    if titrant == "H2SO4":
        titrant_normality = 2
    else:
        titrant_normality = 1
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
            - analyte_mass * (HF + HSO4)
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
    titrant=default.titrant,
    HF=0,
    HSO4=0,
    emf0_guess=None,
):
    """Simple Gran-plot first guesses for alkalinity, EMF0 and pH.
    Uses a subset of the provided data points.
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
    )[0]
    if titrant == "H2SO4":
        alkalinity_guess *= 2
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
                titrant=titrant,
                HF=HF,
                HSO4=HSO4,
            )
        )
    pH_guesses = convert.emf_to_pH(emf, emf0_guess, temperature)
    return alkalinity_guess, emf0_guess, pH_guesses, G


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
    """Calculate residuals for the complete-calculation solver."""
    alkalinity, emf0 = alkalinity_emf0
    pH = convert.emf_to_pH(emf, emf0, temperature)
    mixture_mass = titrant_mass + analyte_mass
    dilution_factor = convert.get_dilution_factor(titrant_mass, analyte_mass)
    residual = (
        simulate.alkalinity(pH, totals, k_constants)
        - alkalinity * dilution_factor
        + titrant_mass * titrant_molinity / mixture_mass
    )
    return residual


def solve_emf_complete(
    titrant_molinity,
    titrant_mass,
    emf,
    temperature,
    analyte_mass,
    totals,
    k_constants,
    least_squares_kwargs=default.least_squares_kwargs,
    pH_range=default.pH_range,
    emf0_guess=None,
):
    """Solve for alkalinity and EMF0 using the complete-calculation method."""
    # Get initial guesses
    alkalinity_guess, emf0_guess, pH_guesses = gran_guesses(
        titrant_mass,
        emf,
        temperature,
        analyte_mass,
        titrant_molinity,
        emf0_guess=emf0_guess,
    )[:-1]
    # Set which data points to use in the final solver
    G = (pH_guesses >= pH_range[0]) & (pH_guesses <= pH_range[1])
    totals_G = {k: v[G] if np.size(v) > 1 else v for k, v in totals.items()}
    k_constants_G = {k: v[G] if np.size(v) > 1 else v for k, v in k_constants.items()}
    # Solve for alkalinity and EMF0
    opt_result = least_squares(
        _lsqfun_solve_emf_complete,
        [alkalinity_guess, emf0_guess],
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
    # Add which data points were used and initial guesses to the output
    opt_result["data_used"] = G
    opt_result["alkalinity_gran"] = alkalinity_guess * 1e6
    opt_result["emf0_gran"] = emf0_guess
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
    """Calculate residuals for the complete-calculation solver."""
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
    least_squares_kwargs=default.least_squares_kwargs,
    pH_range=default.pH_range,
    emf0_guess=None,
):
    """Solve for alkalinity and EMF0 using the complete-calculation method."""
    # Get initial guesses
    alkalinity_guess, emf0_guess, pH_guesses = gran_guesses(
        titrant_mass,
        emf,
        temperature,
        analyte_mass,
        titrant_molinity,
        titrant="H2SO4",
        emf0_guess=emf0_guess,
    )[:-1]
    # Set which data points to use in the final solver
    G = (pH_guesses >= pH_range[0]) & (pH_guesses <= pH_range[1])
    totals_G = {k: v[G] if np.size(v) > 1 else v for k, v in totals.items()}
    k_constants_G = {k: v[G] if np.size(v) > 1 else v for k, v in k_constants.items()}
    # Solve for alkalinity and EMF0
    opt_result = least_squares(
        _lsqfun_solve_emf_complete_H2SO4,
        [alkalinity_guess, emf0_guess],
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
    opt_result["alkalinity_gran"] = alkalinity_guess * 1e6
    opt_result["emf0_gran"] = emf0_guess
    return opt_result


def _lsqfun_calibrate(
    titrant_molinity,
    alkalinity_certified,
    titrant_mass,
    emf,
    temperature,
    analyte_mass,
    totals,
    k_constants,
    solver_kwargs,
):
    """Calculate residuals for the calibrator."""
    alkalinity = solve_emf_complete(
        titrant_molinity[0],
        titrant_mass,
        emf,
        temperature,
        analyte_mass,
        totals,
        k_constants,
        **solver_kwargs,
    )["x"][0]
    return alkalinity * 1e6 - alkalinity_certified


def calibrate(
    alkalinity_certified,
    titrant_mass,
    emf,
    temperature,
    analyte_mass,
    totals,
    k_constants,
    titrant_molinity_guess=None,
    least_squares_kwargs=default.least_squares_kwargs,
    solver_kwargs={},
):
    """Solve for titrant_molinity given alkalinity_certified."""
    if titrant_molinity_guess is None:
        titrant_molinity_guess = copy.deepcopy(default.titrant_molinity_guess)
    solver_kwargs["least_squares_kwargs"] = least_squares_kwargs
    return least_squares(
        _lsqfun_calibrate,
        [titrant_molinity_guess],
        args=(
            alkalinity_certified,
            titrant_mass,
            emf,
            temperature,
            analyte_mass,
            totals,
            k_constants,
            solver_kwargs,
        ),
        **least_squares_kwargs,
    )


def _lsqfun_calibrate_H2SO4(
    titrant_molinity,
    alkalinity_certified,
    titrant_mass,
    emf,
    temperature,
    analyte_mass,
    analyte_total_sulfate,
    analyte_salinity,
    totals,
    k_constants,
    solver_kwargs,
):
    """Calculate residuals for the calibrator."""
    # Update sulfate dilution by titrant and its consequences
    totals["total_sulfate"] = (
        analyte_total_sulfate * analyte_mass + titrant_molinity * titrant_mass
    ) / (analyte_mass + titrant_mass)
    totals_pyco2 = convert.totals_to_pyco2(totals, analyte_salinity)
    k_constants = interface.get_k_constants(totals_pyco2, temperature)
    # Solve for alkalinity
    alkalinity = solve_emf_complete_H2SO4(
        titrant_molinity[0],
        titrant_mass,
        emf,
        temperature,
        analyte_mass,
        totals,
        k_constants,
        **solver_kwargs,
    )["x"][0]
    # Return residual
    return alkalinity * 1e6 - alkalinity_certified


def calibrate_H2SO4(
    alkalinity_certified,
    titrant_mass,
    emf,
    temperature,
    analyte_mass,
    analyte_total_sulfate,
    analyte_salinity,
    totals,
    k_constants,
    titrant_molinity_guess=None,
    least_squares_kwargs=default.least_squares_kwargs,
    solver_kwargs={},
):
    """Solve for titrant_molinity given alkalinity_certified."""
    if titrant_molinity_guess is None:
        titrant_molinity_guess = copy.deepcopy(default.titrant_molinity_guess)
    solver_kwargs["least_squares_kwargs"] = least_squares_kwargs
    return least_squares(
        _lsqfun_calibrate_H2SO4,
        [titrant_molinity_guess],
        args=(
            alkalinity_certified,
            titrant_mass,
            emf,
            temperature,
            analyte_mass,
            analyte_total_sulfate,
            analyte_salinity,
            totals.copy(),
            k_constants,
            solver_kwargs,
        ),
        **least_squares_kwargs,
    )
