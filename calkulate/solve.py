# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew P. Humphreys  (GNU GPLv3)
"""Solve titration data for total alkalinity."""

import numpy as np
from scipy.stats import linregress
from scipy.optimize import least_squares
import PyCO2SYS as pyco2
from . import constants, convert, simulate


def gran_estimator(mixture_mass, emf, temperature, use_points=None):
    """Simple Gran-plot estimator function following DAA03 eq. 10."""
    if use_points is None:
        G = np.full(np.size(emf), True)
    else:
        G = use_points
    temperature_K = temperature[G] + constants.absolute_zero
    return mixture_mass[G] * np.exp(
        emf[G] * constants.faraday / (constants.ideal_gas * temperature_K)
    )


def gran_guess_alkalinity(titration, metadata, use_points=None):
    """Simple Gran-plot first-guess of alkalinity."""
    if use_points is None:
        G = np.full(np.size(titration["titrant_mass"]), True)
    else:
        G = use_points
    mixture_mass = titration["titrant_mass"] + metadata["analyte_mass"]
    gradient, intercept_y = linregress(
        titration["titrant_mass"][G],
        gran_estimator(
            mixture_mass, titration["emf"], titration["temperature"], use_points=G
        ),
    )[:2]
    intercept_x = -intercept_y / gradient
    alkalinity_guess = (
        intercept_x * metadata["titrant_molinity"] / metadata["analyte_mass"]
    )
    return alkalinity_guess


def gran_guess_emf0(titration, metadata, use_points=None, HF=0, HSO4=0):
    """Simple Gran-plot first-guess of EMF0 following DAA03 eq. 11."""
    if use_points is None:
        G = np.full(np.size(titration["titrant_mass"]), True)
    else:
        G = use_points
    alkalinity_guess = gran_guess_alkalinity(titration, metadata, use_points=G)
    temperature_K = titration["temperature"][G] + constants.absolute_zero
    mixture_mass = titration["titrant_mass"] + metadata["analyte_mass"]
    return titration["emf"][G] - (
        constants.ideal_gas * temperature_K / constants.faraday
    ) * np.log(
        (
            (
                titration["titrant_mass"][G] * metadata["titrant_molinity"]
                - metadata["analyte_mass"] * alkalinity_guess
            )
            - metadata["analyte_mass"] * (HF + HSO4)
        )
        / mixture_mass[G]
    )


def gran_guesses(titration, metadata):
    """Simple Gran plot first guesses for alkalinity and EMF0."""
    mixture_mass = titration["titrant_mass"] + metadata["analyte_mass"]
    estimator = gran_estimator(mixture_mass, titration["emf"], titration["temperature"])
    G = (estimator > 0.1 * np.max(estimator)) & (estimator < 0.9 * np.max(estimator))
    alkalinity_guess = gran_guess_alkalinity(titration, metadata, use_points=G)
    emf0_guess = np.mean(gran_guess_emf0(titration, metadata, use_points=G))
    return alkalinity_guess, emf0_guess, G


def _lsqfun_complete_emf(
    alkalinity_emf0, titration, metadata, titrant_molinity,
):
    alkalinity, emf0 = alkalinity_emf0
    pH = convert.emf_to_pH(titration.emf.values, emf0, titration.temperature.values)
    mixture_mass = titration.titrant_mass.values + metadata["analyte_mass"]
    residual = (
        simulate.alkalinity(pH, titration)
        - alkalinity * titration.dilution_factor.values
        + titration.titrant_mass.values * titrant_molinity / mixture_mass
    )
    return residual


def complete_emf(titration, metadata, titrant_molinity=None, pH_range=(3, 4)):
    """Solve for alkalinity and EMF0 using the complete calculation method."""
    if titrant_molinity is None:
        titrant_molinity = metadata["titrant_molinity"]
    alkalinity_guess, emf0_guess = gran_guesses(titration, metadata)[:2]
    pH_guess = convert.emf_to_pH(titration["emf"], emf0_guess, titration["temperature"])
    assert pH_range[0] < pH_range[1]
    G = (pH_guess > pH_range[0]) & (pH_guess < pH_range[1])
    opt_result = least_squares(
        _lsqfun_complete_emf,
        [alkalinity_guess, emf0_guess],
        args=(titration[G], metadata, titrant_molinity),
        method="lm",
        x_scale=[1e-6, 1],
    )
    opt_result["use_points"] = G
    return opt_result


def complete_pH(titration, metadata, titrant_molinity=None, pH_range=(3, 4)):
    """Solve for alkalinity from pH using the complete calculation method."""
    assert pH_range[0] < pH_range[1]
    G = (titration.pH > pH_range[0]) & (titration.pH < pH_range[1])
    if titrant_molinity is None:
        titrant_molinity = metadata["titrant_molinity"]
    alkalinity_points = simulate.alkalinity(titration.pH[G], titration[G])
    submixture_mass = titration["titrant_mass"][G] + metadata["analyte_mass"]
    alkalinity_points += (
        titration["titrant_mass"][G] * titrant_molinity / submixture_mass
    )
    alkalinity_points /= titration.dilution_factor[G]
    return {
        "alkalinity_points": alkalinity_points,
        "x": np.array([np.mean(alkalinity_points)]),
        "alkalinity_std": np.std(alkalinity_points),
        "use_points": G,
    }


def _lsqfun_calibrate(
    titrant_molinity, titration, metadata, solver=complete_emf, pH_range=(3, 4)
):
    alkalinity = (
        solver(
            titration, metadata, titrant_molinity=titrant_molinity, pH_range=pH_range
        )["x"][0]
        * 1e6
    )
    return alkalinity - metadata["alkalinity_certified"]


def calibrate(
    titration,
    metadata,
    titrant_molinity_guess=0.1,
    solver=complete_emf,
    pH_range=(3, 4),
):
    """Calibrate the acid concentration where alkalinity is known."""
    return least_squares(
        _lsqfun_calibrate,
        titrant_molinity_guess,
        args=(titration, metadata),
        kwargs={"solver": solver, "pH_range": pH_range},
        method="lm",
    )
