# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew P. Humphreys  (GNU GPLv3)
"""Solve titration data for total alkalinity."""

import warnings
import numpy as np, PyCO2SYS as pyco2
from scipy.stats import linregress
from scipy.optimize import least_squares
from . import constants, convert, options, simulate

pyco2.solve.get.TAfromTCpH = pyco2.solve.get.TAfromTCpH_fixed


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


def gran_guess_alkalinity(
    titrant_mass, emf, temperature, analyte_mass, titrant_molinity, use_points=None
):
    """Simple Gran-plot first-guess of alkalinity."""
    if use_points is None:
        G = np.full(np.size(titrant_mass), True)
    else:
        G = use_points
    mixture_mass = titrant_mass + analyte_mass
    gradient, intercept_y = linregress(
        titrant_mass[G], gran_estimator(mixture_mass, emf, temperature, use_points=G),
    )[:2]
    intercept_x = -intercept_y / gradient
    alkalinity_guess = intercept_x * titrant_molinity / analyte_mass
    return alkalinity_guess


def gran_guess_emf0(
    titrant_mass,
    emf,
    temperature,
    analyte_mass,
    titrant_molinity,
    use_points=None,
    HF=0,
    HSO4=0,
):
    """Simple Gran-plot first-guess of EMF0 following DAA03 eq. 11."""
    if use_points is None:
        G = np.full(np.size(titrant_mass), True)
    else:
        G = use_points
    alkalinity_guess = gran_guess_alkalinity(
        titrant_mass, emf, temperature, analyte_mass, titrant_molinity, use_points=G
    )
    temperature_K = temperature[G] + constants.absolute_zero
    mixture_mass = titrant_mass + analyte_mass
    return emf[G] - (constants.ideal_gas * temperature_K / constants.faraday) * np.log(
        (
            (titrant_mass[G] * titrant_molinity - analyte_mass * alkalinity_guess)
            - analyte_mass * (HF + HSO4)
        )
        / mixture_mass[G]
    )


def gran_guesses(titrant_mass, emf, temperature, analyte_mass, titrant_molinity):
    """Simple Gran plot first guesses for alkalinity and EMF0."""
    mixture_mass = titrant_mass + analyte_mass
    estimator = gran_estimator(mixture_mass, emf, temperature)
    G = (estimator > 0.1 * np.max(estimator)) & (estimator < 0.9 * np.max(estimator))
    alkalinity_guess = gran_guess_alkalinity(
        titrant_mass, emf, temperature, analyte_mass, titrant_molinity, use_points=G
    )
    emf0_guess = np.mean(
        gran_guess_emf0(
            titrant_mass, emf, temperature, analyte_mass, titrant_molinity, use_points=G
        )
    )
    return alkalinity_guess, emf0_guess, G


def _lsqfun_complete_emf(
    alkalinity_emf0, titration, analyte_mass, titrant_molinity,
):
    tt = titration
    alkalinity, emf0 = alkalinity_emf0
    pH = convert.emf_to_pH(tt["emf"], emf0, tt["temperature"])
    mixture_mass = tt["titrant_mass"] + analyte_mass
    residual = (
        simulate.alkalinity(pH, tt)
        - alkalinity * tt["dilution_factor"]
        + tt["titrant_mass"] * titrant_molinity / mixture_mass
    )
    return residual


def complete_emf(dataset_row, titrant_molinity=None, pH_range=(3, 4)):
    """Solve for alkalinity and EMF0 using the complete calculation method."""
    dsr = dataset_row
    tt = dataset_row["titration"]
    if titrant_molinity is None:
        titrant_molinity = dsr["titrant_molinity"]
    alkalinity_guess, emf0_guess = gran_guesses(
        tt["titrant_mass"],
        tt["emf"],
        tt["temperature"],
        dsr["analyte_mass"],
        titrant_molinity,
    )[:2]
    pH_guess = convert.emf_to_pH(tt["emf"], emf0_guess, tt["temperature"])
    assert pH_range[0] < pH_range[1]
    G = (pH_guess > pH_range[0]) & (pH_guess < pH_range[1])
    opt_result = least_squares(
        _lsqfun_complete_emf,
        [alkalinity_guess, emf0_guess],
        args=(tt[G], dsr["analyte_mass"], titrant_molinity),
        method=options.solver_method,
        x_scale=[1e-6, 1],
    )
    opt_result["use_points"] = G
    return opt_result


def complete_pH(dataset_row, titrant_molinity=None, pH_range=(3, 4)):
    """Solve for alkalinity from pH using the complete calculation method."""
    dsr = dataset_row
    tt = dataset_row["titration"]
    assert pH_range[0] < pH_range[1]
    G = (tt["pH"] > pH_range[0]) & (tt["pH"] < pH_range[1])
    if titrant_molinity is None:
        titrant_molinity = dsr["titrant_molinity"]
    alkalinity_points = simulate.alkalinity(tt["pH"][G], tt[G])
    submixture_mass = tt["titrant_mass"][G] + dsr["analyte_mass"]
    alkalinity_points += tt["titrant_mass"][G] * titrant_molinity / submixture_mass
    alkalinity_points /= tt["dilution_factor"][G]
    return {
        "alkalinity_points": alkalinity_points,
        "x": np.array([np.mean(alkalinity_points)]),
        "alkalinity_std": np.std(alkalinity_points),
        "use_points": G,
    }


def _get_PyCO2SYS_inputs(titration, salinity):
    totals = {
        convert.calk_to_pyco2[t]: titration[t].values * 1e-6
        for t in [
            "total_silicate",
            "total_phosphate",
            "total_ammonia",
            "total_sulfide",
            "total_sulfate",
            "total_borate",
            "total_fluoride",
        ]
    }
    k_constants = {
        convert.calk_to_pyco2[k]: titration[k].values
        for k in [
            "k_ammonia",
            "k_borate",
            "k_bisulfate",
            "k_carbonic_1",
            "k_carbonic_2",
            "k_fluoride",
            "k_phosphate_1",
            "k_phosphate_2",
            "k_phosphate_3",
            "k_silicate",
            "k_sulfide",
            "k_water",
        ]
    }
    k_constants = pyco2.convert.get_pHfactor_to_Free(
        titration["temperature"].values,
        salinity,
        totals,
        k_constants,
        3,
        options.pyco2_opt_k_carbonic,
    )
    return totals, k_constants


def _lsqfun_calibrate_emf(
    titrant_molinity__emf0,
    titrant_mass,
    emf,
    temperature,
    dic,
    analyte_mass,
    alkalinity_certified,
    totals,
    k_constants,
):
    titrant_molinity, emf0 = titrant_molinity__emf0
    mixture_alkalinity = (
        alkalinity_certified * analyte_mass - titrant_mass * titrant_molinity * 1e6
    ) / (analyte_mass + titrant_mass)
    pH_emf0 = -np.log10(convert.emf_to_h(emf, emf0, temperature))
    pH_alk = pyco2.solve.get.pHfromTATC(
        mixture_alkalinity * 1e-6, dic * 1e-6, totals, k_constants
    )
    return pH_emf0 - pH_alk


def calibrate_emf(dataset_row, titrant_molinity_guess=0.1, pH_range=(3, 4)):
    """Calibrate the titrant molinity where alkalinity is known for EMF data."""
    dsr = dataset_row
    tt = dataset_row["titration"]
    mixture_mass = tt["titrant_mass"] + dsr["analyte_mass"]
    estimator = gran_estimator(mixture_mass, tt["emf"], tt["temperature"])
    G = (estimator > 0.1 * np.max(estimator)) & (estimator < 0.9 * np.max(estimator))
    emf0_guess = np.mean(
        gran_guess_emf0(
            tt["titrant_mass"],
            tt["emf"],
            tt["temperature"],
            dsr["analyte_mass"],
            titrant_molinity_guess,
            use_points=G,
        )
    )
    totals, k_constants = _get_PyCO2SYS_inputs(tt[G], dsr["salinity"])
    return least_squares(
        _lsqfun_calibrate_emf,
        [titrant_molinity_guess, emf0_guess],
        args=(
            tt["titrant_mass"][G].to_numpy(),
            tt["emf"][G].to_numpy(),
            tt["temperature"][G].to_numpy(),
            tt["dic"][G].to_numpy(),
            dsr["analyte_mass"],
            dsr["alkalinity_certified"],
            totals,
            k_constants,
        ),
        method=options.solver_method,
    )


def _lsqfun_calibrate_pH(
    titrant_molinity,
    titrant_mass,
    pH,
    dic,
    analyte_mass,
    alkalinity_certified,
    totals,
    k_constants,
):
    mixture_alkalinity = (
        alkalinity_certified * analyte_mass - titrant_mass * titrant_molinity * 1e6
    ) / (analyte_mass + titrant_mass)
    pH_alk = pyco2.solve.get.pHfromTATC(
        mixture_alkalinity * 1e-6, dic * 1e-6, totals, k_constants
    )
    return pH - pH_alk


def calibrate_pH(dataset_row, titrant_molinity_guess=0.1, pH_range=(3, 4)):
    """Calibrate the titrant molinity where alkalinity is known for pH data."""
    dsr = dataset_row
    tt = dataset_row["titration"]
    assert pH_range[0] < pH_range[1]
    G = (tt["pH"] > pH_range[0]) & (tt["pH"] < pH_range[1])
    assert sum(G) > 0, "No titration points fall in the specified pH range!"
    totals, k_constants = _get_PyCO2SYS_inputs(tt[G], dsr["salinity"])
    return least_squares(
        _lsqfun_calibrate_pH,
        titrant_molinity_guess,
        args=(
            tt["titrant_mass"][G].to_numpy(),
            tt["pH"][G].to_numpy(),
            tt["dic"][G].to_numpy(),
            dsr["analyte_mass"],
            dsr["alkalinity_certified"],
            totals,
            k_constants,
        ),
        method=options.solver_method,
    )


def _lsqfun_calibrate(titrant_molinity, dataset_row, solver, pH_range):
    alkalinity = solver(
        dataset_row, titrant_molinity=titrant_molinity, pH_range=pH_range
    )["x"][0]
    return alkalinity * 1e6 - dataset_row["alkalinity_certified"]


def calibrate(
    dataset_row, titrant_molinity_guess=0.1, pH_range=(3, 4), solver=complete_emf
):
    """Calibrate the titrant molinity where alkalinity is known."""
    warnings.filterwarnings("ignore")
    return least_squares(
        _lsqfun_calibrate,
        titrant_molinity_guess,
        args=(dataset_row, solver, pH_range),
        method=options.solver_method,
    )
