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


def gran_guess_alkalinity(titrant, analyte, emf, temperature, use_points=None):
    """Simple Gran-plot first-guess of alkalinity."""
    if use_points is None:
        G = np.full(np.size(emf), True)
    else:
        G = use_points
    mixture_mass = titrant["mass"] + analyte["mass"]
    gradient, intercept_y = linregress(
        titrant["mass"][G], gran_estimator(mixture_mass, emf, temperature, use_points=G)
    )[:2]
    intercept_x = -intercept_y / gradient
    alkalinity_guess = intercept_x * titrant["molinity"] / analyte["mass"]
    return alkalinity_guess


def gran_guess_emf0(titrant, analyte, emf, temperature, use_points=None, HF=0, HSO4=0):
    """Simple Gran-plot first-guess of EMF0 following DAA03 eq. 11."""
    if use_points is None:
        G = np.full(np.size(emf), True)
    else:
        G = use_points
    alkalinity_guess = gran_guess_alkalinity(
        titrant, analyte, emf, temperature, use_points=G
    )
    temperature_K = temperature[G] + constants.absolute_zero
    mixture_mass = titrant["mass"] + analyte["mass"]
    return emf[G] - (constants.ideal_gas * temperature_K / constants.faraday) * np.log(
        (
            (
                titrant["mass"][G] * titrant["molinity"]
                - analyte["mass"] * alkalinity_guess
            )
            - analyte["mass"] * (HF + HSO4)
        )
        / mixture_mass[G]
    )


def gran_guesses(titrant, analyte, emf, temperature):
    """Simple Gran plot first guesses for alkalinity and EMF0."""
    mixture_mass = titrant["mass"] + analyte["mass"]
    estimator = gran_estimator(mixture_mass, emf, temperature)
    G = (estimator > 0.1 * np.max(estimator)) & (estimator < 0.9 * np.max(estimator))
    alkalinity_guess = gran_guess_alkalinity(
        titrant, analyte, emf, temperature, use_points=G
    )
    emf0_guess = np.mean(
        gran_guess_emf0(titrant, analyte, emf, temperature, use_points=G)
    )
    return alkalinity_guess, emf0_guess, G


def _lsqfun_complete_emf(
    alkalinity_emf0,
    titrant,
    analyte,
    emf,
    temperature,
    totals,
    k_constants,
    dilution_factor=1,
):
    alkalinity, emf0 = alkalinity_emf0
    pH = convert.emf_to_pH(emf, emf0, temperature)
    mixture_mass = titrant["mass"] + analyte["mass"]
    return (
        simulate.alkalinity(pH, totals, k_constants, dic=analyte["dic"] * 1e-6)
        - alkalinity * dilution_factor
        + titrant["mass"] * titrant["molinity"] / mixture_mass
    )


def subset(d, use_points):
    return {
        k: v[use_points] if np.shape(v) == np.shape(use_points) else v
        for k, v in d.items()
    }


def complete_emf(
    titrant, analyte, emf, temperature, totals, k_constants, pH_range=(3, 4)
):
    """Solve for alkalinity and EMF0 using the complete calculation method."""
    alkalinity_guess, emf0_guess = gran_guesses(titrant, analyte, emf, temperature)[:2]
    pH_guess = convert.emf_to_pH(emf, emf0_guess, temperature)
    assert pH_range[0] < pH_range[1]
    G = (pH_guess > pH_range[0]) & (pH_guess < pH_range[1])
    opt_result = least_squares(
        _lsqfun_complete_emf,
        [alkalinity_guess, emf0_guess],
        args=(
            subset(titrant, G),
            subset(analyte, G),
            emf[G],
            temperature[G],
            subset(totals, G),
            subset(k_constants, G),
        ),
        kwargs={
            "dilution_factor": convert.dilution_factor(
                analyte["mass"], analyte["mass"] + titrant["mass"][G]
            ),
        },
        method="lm",
        x_scale=[1e-6, 1],
    )
    opt_result["use_points"] = G
    return opt_result


def complete_pH(
    titrant, analyte, pH, temperature, totals, k_constants, pH_range=(3, 4)
):
    """Solve for alkalinity from pH using the complete calculation method."""
    assert pH_range[0] < pH_range[1]
    G = (pH > pH_range[0]) & (pH < pH_range[1])
    alkalinity_points = simulate.alkalinity(
        pH, totals, k_constants, dic=analyte["dic"] * 1e-6
    )[G]
    submixture_mass = titrant["mass"][G] + analyte["mass"]
    alkalinity_points += titrant["mass"][G] * titrant["molinity"] / submixture_mass
    alkalinity_points /= convert.dilution_factor(analyte["mass"], submixture_mass)
    return {
        "alkalinity_points": alkalinity_points,
        "x": np.array([np.mean(alkalinity_points)]),
        "alkalinity_std": np.std(alkalinity_points),
        "use_points": G,
    }


def _lsqfun_calibrate(
    titrant_molinity,
    titrant,
    analyte,
    emf_or_pH,
    temperature,
    totals,
    k_constants,
    solver=complete_emf,
    pH_range=(3, 4),
):
    titrant["molinity"] = titrant_molinity
    alkalinity = (
        solver(
            titrant,
            analyte,
            emf_or_pH,
            temperature,
            totals,
            k_constants,
            pH_range=pH_range,
        )["x"][0]
        * 1e6
    )
    return alkalinity - analyte["alkalinity_certified"]


def calibrate(
    titrant,
    analyte,
    emf_or_pH,
    temperature,
    totals,
    k_constants,
    titrant_molinity_guess=0.1,
    solver=complete_emf,
    pH_range=(3, 4),
):
    """Calibrate the acid concentration where alkalinity is known."""
    return least_squares(
        _lsqfun_calibrate,
        titrant_molinity_guess,
        args=(titrant, analyte, emf_or_pH, temperature, totals, k_constants),
        kwargs={"pH_range": pH_range},
        method="lm",
    )
