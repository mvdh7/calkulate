# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew P. Humphreys  (GNU GPLv3)
"""Solve titration data for total alkalinity."""

import numpy as np
from scipy.stats import linregress
from scipy.optimize import least_squares
import PyCO2SYS as pyco2
from . import constants, convert


def gran_estimator(titration, use_points=None):
    """Simple Gran-plot estimator function following DAA03 eq. 10."""
    tt = titration
    G = use_points
    if G is None:
        G = np.full(np.size(titration.mixture.temperature), True)
    temperature_K = tt.mixture.temperature[G] + constants.absolute_zero
    return tt.mixture.mass[G] * np.exp(
        tt.mixture.emf[G] * constants.faraday / (constants.ideal_gas * temperature_K)
    )


def gran_guess_alkalinity(titration, use_points=None):
    """Simple Gran-plot first-guess of alkalinity."""
    tt = titration
    G = use_points
    if G is None:
        G = np.full(np.size(titration.mixture.temperature), True)
    gradient, intercept_y, _, _, _ = linregress(
        tt.titrant.mass[G], gran_estimator(tt, use_points=G)
    )
    intercept_x = -intercept_y / gradient
    alkalinity_guess = intercept_x * tt.titrant.concentration / tt.analyte.mass
    return alkalinity_guess


def gran_guess_emf0(titration, use_points=None, HF=0, HSO4=0):
    """Simple Gran-plot first-guess of EMF0 following DAA03 eq. 11."""
    tt = titration
    G = use_points
    if G is None:
        G = np.full(np.size(titration.mixture.temperature), True)
    alkalinity_guess = gran_guess_alkalinity(tt, use_points=G)
    temperature_K = tt.mixture.temperature[G] + constants.absolute_zero
    return tt.mixture.emf[G] - (
        constants.ideal_gas * temperature_K / constants.faraday
    ) * np.log(
        (
            (
                tt.titrant.mass[G] * tt.titrant.molinity
                - tt.analyte.mass * alkalinity_guess
            )
            - tt.analyte.mass * (HF + HSO4)
        )
        / tt.mixture.mass[G]
    )


def gran_guesses(titration):
    """Simple Gran plot first guesses for alkalinity and EMF0."""
    tt = titration
    estimator = gran_estimator(tt)
    G = (estimator > 0.1 * np.max(estimator)) & (estimator < 0.9 * np.max(estimator))
    alkalinity_guess = gran_guess_alkalinity(tt, use_points=G)
    emf0_guess = np.mean(gran_guess_emf0(tt, use_points=G))
    return alkalinity_guess, emf0_guess, G


def _lsqfun_complete(alkalinity_emf0, mixture, titrant):
    alkalinity, emf0 = alkalinity_emf0
    pH = convert.emf_to_pH(mixture.emf, emf0, mixture.temperature)
    return (
        pyco2.solve.get.TAfromTCpH(
            mixture.total_carbonate * 1e-6,
            pH,
            mixture.total_salts,
            mixture.equilibrium_constants,
        )
        - alkalinity * mixture.dilution_factor
        + titrant.mass * titrant.molinity / mixture.mass
    )


def complete(titration):
    """Solve for alkalinity and EMF0 using the complete calculation method."""
    tt = titration
    pH_guess = convert.emf_to_pH(
        tt.mixture.emf, tt.analyte.emf0_guess, tt.mixture.temperature
    )
    use_points = (pH_guess > 3) & (pH_guess < 4)
    opt_result = least_squares(
        _lsqfun_complete,
        [tt.analyte.alkalinity_guess * 1e-6, tt.analyte.emf0_guess],
        args=(tt.mixture.subset(use_points), tt.titrant.subset(use_points)),
        method="lm",
        x_scale=[1e-6, 1],
    )
    opt_result["use_points"] = use_points
    return opt_result
