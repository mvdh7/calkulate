# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019--2022  Matthew P. Humphreys  (GNU GPLv3)
"""Simulate solution properties during a titration."""

import numpy as np
import PyCO2SYS as pyco2
from . import convert, default, io
from .titration import Titration


def alkalinity_components(pH, totals, k_constants, opt_pH_scale=default.opt_pH_scale):
    """Calculate chemical speciation from pH.

    totals should include dilution correction and be in mol/kg-solution.

    k_constants should be on the scale specified by opt_pH_scale, with:
        1 = Total                pH = -log10([H+] + [HSO4-])
        2 = Seawater             pH = -log10([H+] + [HSO4-] + [HF])
        3 = Free      [default]  pH = -log10([H+])
    Note that when using the Total scale, k_fluoride must also be on the Total scale
    (this does not happen by default in PyCO2SYS!).

    Outputs are substance contents in mol/kg-solution.
    """
    # Check opt_pH_scale is valid
    opt_pH_scales = [1, 2, 3]
    assert (
        opt_pH_scale in opt_pH_scales
    ), "opt_pH_scale must be 1 (Total), 2 (Seawater) or 3 (Free)."
    # Build up dict of solution components
    components = {}
    h = components["H"] = 10.0 ** -pH
    if "k_water" in k_constants:
        components["OH"] = k_constants["k_water"] / h
    if "dic" in totals:
        TCO2 = totals["dic"]
        K1 = k_constants["k_carbonic_1"]
        K2 = k_constants["k_carbonic_2"]
        components["CO2"] = TCO2 / (1 + K1 / h + K1 * K2 / h ** 2)
        components["HCO3"] = K1 * components["CO2"] / h
        components["CO3"] = K2 * components["HCO3"] / h
    if "total_borate" in totals:
        TB = totals["total_borate"]
        KB = k_constants["k_borate"]
        components["BOH4"] = TB * KB / (KB + h)
    if "total_phosphate" in totals:
        TPO4 = totals["total_phosphate"]
        KP1 = k_constants["k_phosphoric_1"]
        KP2 = k_constants["k_phosphoric_2"]
        KP3 = k_constants["k_phosphoric_3"]
        phosphoric_denom = h ** 3 + KP1 * h ** 2 + KP1 * KP2 * h + KP1 * KP2 * KP3
        components["H3PO4"] = TPO4 * h ** 3 / phosphoric_denom
        components["HPO4"] = TPO4 * KP1 * KP2 * h / phosphoric_denom
        components["PO4"] = TPO4 * KP1 * KP2 * KP3 / phosphoric_denom
    if "total_silicate" in totals:
        TSi = totals["total_silicate"]
        KSi = k_constants["k_silicate"]
        components["H3SiO4"] = TSi * KSi / (KSi + h)
    if "total_ammonia" in totals:
        TNH3 = totals["total_ammonia"]
        KNH3 = k_constants["k_ammonia"]
        components["NH3"] = TNH3 * KNH3 / (KNH3 + h)
    if "total_sulfide" in totals:
        TH2S = totals["total_sulfide"]
        KH2S = k_constants["k_sulfide"]
        components["HS"] = TH2S * KH2S / (KH2S + h)
    if "total_alpha" in totals:
        total_alpha = totals["total_alpha"]
        k_alpha = k_constants["k_alpha"]
        alpha = total_alpha * k_alpha / (k_alpha + h)
        alphaH = total_alpha - alpha
        components["alk_alpha"] = np.where(
            -np.log10(k_alpha) <= default.zlp, -alphaH, alpha
        )
    if "total_beta" in totals:
        total_beta = totals["total_beta"]
        k_beta = k_constants["k_beta"]
        beta = total_beta * k_beta / (k_beta + h)
        betaH = total_beta - beta
        components["alk_beta"] = np.where(
            -np.log10(k_beta) <= default.zlp, -betaH, beta
        )
    # pH-scale-dependent components
    if opt_pH_scale in [1, 3]:
        if "total_fluoride" in totals:
            TF = totals["total_fluoride"]
            KF = k_constants["k_fluoride"]
            components["HF"] = TF * h / (KF + h)
    if opt_pH_scale == 3:
        if "total_sulfate" in totals:
            TSO4 = totals["total_sulfate"]
            KSO4 = k_constants["k_bisulfate"]
            components["HSO4"] = TSO4 * h / (KSO4 + h)
    return components


# Multipliers for each component in the alkalinity equation.
# Keys correspond to the output dict of alkalinity_components().
component_multipliers = {
    "H": -1,
    "OH": +1,
    "CO2": 0,
    "HCO3": +1,
    "CO3": +2,
    "BOH4": +1,
    "HSO4": -1,
    "HF": -1,
    "H3PO4": -1,
    "HPO4": +1,
    "PO4": +2,
    "H3SiO4": +1,
    "NH3": +1,
    "HS": +1,
    "alk_alpha": +1,
    "alk_beta": +1,
    "alkalinity_estimate": +1,  # for plotting
}


def alkalinity(pH, totals, k_constants, opt_pH_scale=default.opt_pH_scale):
    """Estimate total alkalinity from pH and total salts in mol/kg-solution."""
    components = alkalinity_components(
        pH, totals, k_constants, opt_pH_scale=opt_pH_scale
    )
    return np.sum([v * component_multipliers[k] for k, v in components.items()], axis=0)


def _titration(
    alkalinity,
    analyte_mass=0.1,
    dic=0,
    emf0=600,
    salinity=35,
    temperature=25,
    titrant_mass_start=0,
    titrant_mass_step=0.15e-3,
    titrant_mass_stop=4.2e-3,
    titrant_molinity=0.1,
    **pyco2sys_kwargs,
):
    """Simulate an titration of a seawater analyte with HCl.

    Parameters
    ----------
    alkalinity : float
        Total alkalinity content of the analyte in µmol/kg.
    analyte_mass : float, optional
        Mass of the analyte in kg, by default 0.1 kg.
    dic : float, optional
        Dissolved inorganic carbon of the analyte in µmol/kg, by default 0 µmol/kg.
    emf0 : float, optional
        EMF0 of the electrode in mV, by default 600 mV.
    salinity : float, optional
        Practical salinity of the analyte, by default 35.
    temperature : float, optional
        Temperature of the analyte in °C, by default 25 °C.
    titrant_mass_start : float, optional
        Mass of titrant at the start of the titration in kg, by default 0 kg.
    titrant_mass_step : float, optional
        Mass of each titrant addition step in kg, by default 0.15e-3 kg.
    titrant_mass_stop : float, optional
        Mass at which to stop the titration (exclusive) in kg, by default 4.2e-3 kg.
    titrant_molinity : float, optional
        Molinity of the titrant in mol/kg, by default 0.1 mol/kg.
    **pyco2sys_kwargs
        Additional kwargs passed on to PyCO2SYS.

    Returns
    -------
    titrant_mass : array_like
        Mass of the titrant through the titration in kg.
    emf : array_like
        EMF across the titrant-analyte mixture through the titration in mV.
    temperature : array_like
        Temperature through the titration in °C.
    analyte_mass : float
        Mass of the analyte in kg.
    totals : dict of array_like
        Total salt contents through the titration in mol/kg.
    k_constants : dict of array_like
        Stoichiometric equilibrium constants through the titration.
    """
    # Create arrays of titrant_mass in kg and temperature in °C
    titrant_mass = np.arange(titrant_mass_start, titrant_mass_stop, titrant_mass_step)
    if np.isscalar(temperature):
        temperature = np.full_like(titrant_mass, temperature)
    # Ensure we use Calkulate's default PyCO2SYS options if they're not provided
    # within pyco2sys_kwargs
    for opt in [
        "opt_gas_constant",
        "opt_k_bisulfate",
        "opt_k_carbonic",
        "opt_k_fluoride",
        "opt_total_borate",
    ]:
        if opt not in pyco2sys_kwargs:
            pyco2sys_kwargs[opt] = getattr(default, opt)
    # Set up dict of start-condition kwargs and then get totals from PyCO2SYS
    kwargs_start = dict(
        salinity=salinity,
        temperature=temperature,
        opt_pH_scale=3,
        **pyco2sys_kwargs,
    )
    co2sys_start = pyco2.sys(**kwargs_start)
    # Calculate dilution through the titration, keeping units in µmol/kg
    dilution_factor = analyte_mass / (analyte_mass + titrant_mass)
    alkalinity_titration = (
        1e6
        * (analyte_mass * alkalinity * 1e-6 - titrant_mass * titrant_molinity)
        / (analyte_mass + titrant_mass)
    )
    dic_titration = dic * dilution_factor
    kwargs_titration = kwargs_start.copy()
    for k in [
        "total_alpha",
        "total_ammonia",
        "total_beta",
        "total_borate",
        "total_calcium",
        "total_fluoride",
        "total_phosphate",
        "total_silicate",
        "total_sulfate",
        "total_sulfide",
    ]:
        kwargs_titration[k] = co2sys_start[k] * dilution_factor
    # Simulate the titration pH and convert it to EMF
    co2sys_titration = pyco2.sys(
        alkalinity_titration, dic_titration, 1, 2, **kwargs_titration
    )
    pH_titration = co2sys_titration["pH_free"]
    emf = convert.pH_to_emf(pH_titration, emf0, temperature)
    # Get totals (in mol/kg) and k_constants dicts for other Calkulate functions
    totals = {
        k: v * 1e-6 for k, v in co2sys_titration.items() if k.startswith("total_")
    }
    totals["dic"] = co2sys_titration["dic"] * 1e-6
    k_constants = {k: v for k, v in co2sys_titration.items() if k.startswith("k_")}
    return titrant_mass, emf, temperature, analyte_mass, totals, k_constants


def titration(
    alkalinity,
    analyte_mass=0.1,
    dic=0,
    emf0=600,
    salinity=35,
    temperature=25,
    titrant_mass_start=0,
    titrant_mass_step=0.15e-3,
    titrant_mass_stop=4.2e-3,
    titrant_molinity=0.1,
    least_squares_kwargs=default.least_squares_kwargs,
    pH_range=default.pH_range,
    **pyco2sys_kwargs,
):
    """Simulate a titration and return a calibrated and solved Titration object.

    Parameters
    ----------
    alkalinity : float
        Total alkalinity content of the analyte in µmol/kg.
    analyte_mass : float, optional
        Mass of the analyte in kg, by default 0.1 kg.
    dic : float, optional
        Dissolved inorganic carbon of the analyte in µmol/kg, by default 0 µmol/kg.
    emf0 : float, optional
        EMF0 of the electrode in mV, by default 600 mV.
    salinity : float, optional
        Practical salinity of the analyte, by default 35.
    temperature : float, optional
        Temperature of the analyte in °C, by default 25 °C.
    titrant_mass_start : float, optional
        Mass of titrant at the start of the titration in kg, by default 0 kg.
    titrant_mass_step : float, optional
        Mass of each titrant addition step in kg, by default 0.15e-3 kg.
    titrant_mass_stop : float, optional
        Mass at which to stop the titration (exclusive) in kg, by default 4.2e-3 kg.
    titrant_molinity : float, optional
        Molinity of the titrant in mol/kg, by default 0.1 mol/kg.
    least_squares_kwargs : dict, optional
        Additional kwargs passed on to the least-squares solver.
    pH_range : tuple, optional
        Range of pH values to determine alkalinity within, by default (3, 4).
    **pyco2sys_kwargs
        Additional kwargs passed on to PyCO2SYS.

    Returns
    -------
    calkulate.Titration
        A self-calibrated and solved titration dataset.
    """
    tt = Titration(
        salinity=salinity,
        analyte_mass=analyte_mass,
        simulate_alkalinity=alkalinity,
        simulate_kwargs=dict(
            dic=dic,
            emf0=emf0,
            salinity=salinity,
            temperature=temperature,
            titrant_mass_start=titrant_mass_start,
            titrant_mass_step=titrant_mass_step,
            titrant_mass_stop=titrant_mass_stop,
            titrant_molinity=titrant_molinity,
            **pyco2sys_kwargs,
        ),
    )
    tt.calkulate(
        alkalinity,
        analyte_total_sulfate=None,  # H2SO4 simulations not implemented yet
        least_squares_kwargs=least_squares_kwargs,
        pH_range=pH_range,
        titrant_molinity_guess=titrant_molinity,
        titrant="HCl",  # H2SO4 simulations not implemented yet
    )
    return tt
