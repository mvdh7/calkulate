# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew P. Humphreys  (GNU GPLv3)
"""Solve titration data for total alkalinity."""

import numpy as np
from scipy.stats import linregress
from scipy.optimize import least_squares
import PyCO2SYS as pyco2
from . import constants, convert, simulate


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
    alkalinity_guess = intercept_x * tt.titrant.molinity / tt.analyte.mass
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


def _lsqfun_complete_emf(alkalinity_emf0, mixture, titrant):
    alkalinity, emf0 = alkalinity_emf0
    pH = convert.emf_to_pH(mixture.emf, emf0, mixture.temperature)
    return (
        simulate.alkalinity(
            pH,
            mixture.total_salts,
            mixture.equilibrium_constants,
            total_carbonate=mixture.total_carbonate * 1e-6,
        )
        - alkalinity * mixture.dilution_factor
        + titrant.mass * titrant.molinity / mixture.mass
    )


def complete_emf(titration, pH_range=(3, 4)):
    """Solve for alkalinity and EMF0 using the complete calculation method."""
    tt = titration
    alkalinity_guess, tt.analyte.emf0_guess, _ = gran_guesses(tt)
    tt.analyte.alkalinity_guess = alkalinity_guess * 1e6
    pH_guess = convert.emf_to_pH(
        tt.mixture.emf, tt.analyte.emf0_guess, tt.mixture.temperature
    )
    assert pH_range[0] < pH_range[1]
    use_points = (pH_guess > pH_range[0]) & (pH_guess < pH_range[1])
    opt_result = least_squares(
        _lsqfun_complete_emf,
        [tt.analyte.alkalinity_guess * 1e-6, tt.analyte.emf0_guess],
        args=(tt.mixture.subset(use_points), tt.titrant.subset(use_points)),
        method="lm",
        x_scale=[1e-6, 1],
    )
    opt_result["use_points"] = use_points
    return opt_result


def complete_pH(titration, pH_range=(3, 4)):
    """Solve for alkalinity from pH using the complete calculation method."""
    tt = titration
    assert pH_range[0] < pH_range[1]
    use_points = (tt.mixture.pH > pH_range[0]) & (tt.mixture.pH < pH_range[1])
    submixture = tt.mixture.subset(use_points=use_points)
    subtitrant = tt.titrant.subset(use_points=use_points)
    # # v1: PyCO2SYS - is there a problem with its alkalinity equation?!
    # alkalinity_points = pyco2.solve.get.TAfromTCpH(
    #     submixture.total_carbonate * 1e-6,
    #     submixture.pH,
    #     submixture.total_salts,
    #     submixture.equilibrium_constants,
    # )
    # free_to_total = pyco2.convert.free2tot(
    #     submixture.total_salts["TSO4"], submixture.equilibrium_constants["KSO4"]
    # )
    # alkalinity_components = pyco2.solve.get.AlkParts(
    #     submixture.total_carbonate * 1e-6,
    #     submixture.pH,
    #     free_to_total,
    #     submixture.total_salts,
    #     submixture.equilibrium_constants,
    # )
    # h = 10.0 ** -submixture.pH
    # alkalinity_points += alkalinity_components["Hfree"] - h
    # hso4 = (
    #     submixture.total_salts["TSO4"]
    #     * h
    #     / (submixture.equilibrium_constants["KSO4"] + h)
    # )
    # alkalinity_points += alkalinity_components["HSO4"] - hso4
    # hf = submixture.total_salts["TF"] * h / (submixture.equilibrium_constants["KF"] + h)
    # alkalinity_points += alkalinity_components["HF"] - hf
    # v2: calk.simulate
    alkalinity_points = simulate.alkalinity(
        submixture.pH,
        submixture.total_salts,
        submixture.equilibrium_constants,
        total_carbonate=submixture.total_carbonate * 1e-6,
    )
    alkalinity_points += subtitrant.mass * subtitrant.molinity / submixture.mass
    alkalinity_points /= submixture.dilution_factor
    return {
        "alkalinity_points": alkalinity_points,
        "x": np.array([np.mean(alkalinity_points)]),
        "alkalinity_std": np.std(alkalinity_points),
        "use_points": use_points,
    }


def _lsqfun_calibrate(titrant_molinity, titration, solver):
    tt = titration
    tt.titrant.molinity = titrant_molinity[0]
    alkalinity = solver(titration)["x"][0] * 1e6
    return alkalinity - tt.analyte.alkalinity_certified


def calibrate(titration, solver=complete_emf, x0=0.1):
    """Calibrate the acid concentration where alkalinity is known."""
    return least_squares(_lsqfun_calibrate, x0, args=(titration, solver), method="lm")
