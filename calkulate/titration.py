# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019--2025  Matthew P. Humphreys  (GNU GPLv3)
"""Work with titration data in a file.

This is the middle level of the three layers of processing functions.
It includes convenience functions for importing data files and reformatting
them to use with the `core` functions, one file at a time.

The processing steps are
  1.  `read.read_dat`: import titration data file.
  2.  `convert_amount_units`: convert amount units to kg (titrant and analyte).
  3.  `get_totals_k_constants`: calculate total salt contents and equilibrium
      constants.
  4a. `calibrate`: find best fitting titrant molinity given a certified
      alkalinity value.
  4b. `solve`: find best fitting alkalinity and EMF0 given titrant molinity.
"""

import inspect
from collections import namedtuple

import numpy as np
import pandas as pd

from . import convert, core, density, interface


Converted = namedtuple(
    "ConvertResult",
    (
        "titrant_mass",
        "measurement",
        "temperature",
        "analyte_mass",
    ),
)
SolveResult = namedtuple(
    "SolveResult",
    (
        "alkalinity",
        "emf0",
        "pH_initial",
        "temperature_initial",
        "analyte_mass",
        "opt_result",
    ),
)


def convert_amount_units(
    dat_data,
    salinity,
    analyte_mass=None,
    analyte_volume=None,
    molinity_H2SO4=0.1,
    molinity_HCl=0.1,
    molinity_NaCl=0.6,
    temperature_override=None,
    titrant_amount_unit="ml",
    titrant_density=None,
    titrant="HCl",
):
    """Convert titrant and analyte units to mass in kg.

    Parameters
    ----------
    dat_data : DatData
        Data from a titration file imported with `read.read_dat`.
    salinity : float
        Practical salinity of the analyte.
    analyte_mass : _type_, optional
        Analyte mass in kg.  Either this or `analyte_volume` must be given.
    analyte_volume : _type_, optional
        Analyte volume in ml, which is converted to kg assuming analyte is
        seawater.  Either this or `analyte_mass` must be given.
    molinity_H2SO4 : float, optional
        H2SO4 titrant molinity in mol/kg-solution for density calculation, by
        default 0.1.
    molinity_HCl : float, optional
        HCl titrant molinity in mol/kg-solution for density calculation, by
        default 0.1.
    molinity_NaCl : float, optional
        NaCl molinity in mol/kg-solution in the titrant where it is an HCl-NaCl
        mixture, for density calculation, by default 0.6.
    temperature_override : float, optional
        A temperature in °C to use instead of the temperature data in the
        titration file, by default `None`.
    titrant_amount_unit : str, optional
        Units for the titrant amount in the file, one of "ml" (default), "g" or
        "kg".
    titrant_density : float, optional
        Titrant density in g/ml, by default `None`, in which case it is
        calculated from the molinities provided.
    titrant : str, optional
        Which titrant was used, "HCl" (default) or "H2SO4".

    Returns
    -------
    Converted - namedtuple containing the fields
        titrant_mass : array-like
            Titrant mass through the titration in kg.
        measurement : array-like
            EMF through the titration in mV, or pH.
        temperature : array-like
            Temperature through the titration in °C.
        analyte_mass : float
            Mass of the analyte in kg.
    """
    dd = dat_data
    # Overwrite temperature, if requested
    if temperature_override is not None:
        temperature = np.full_like(dd.temperature, temperature_override)
    else:
        temperature = dd.temperature
    # Get titrant mass
    assert titrant_amount_unit.lower() in ["ml", "g", "kg"]
    if titrant_amount_unit.lower() == "ml":
        if titrant_density is None:
            assert titrant.upper() in ["H2SO4", "HCL"]
            if titrant.upper() == "H2SO4":
                titrant_mass = (
                    dd.titrant_amount
                    * density.H2SO4_25C_EAIM(molinity_H2SO4)
                    * 1e-3
                )
            else:
                titrant_mass = (
                    dd.titrant_amount
                    * density.HCl_NaCl_25C_DSC07(
                        molinity_HCl=molinity_HCl,
                        molinity_NaCl=molinity_NaCl,
                    )
                    * 1e-3
                )
        else:
            titrant_mass = dd.titrant_amount * titrant_density * 1e-3
    elif titrant_amount_unit.lower() == "g":
        titrant_mass = dd.titrant_amount * 1e-3
    elif titrant_amount_unit.lower() == "kg":
        titrant_mass = dd.titrant_amount
    # Convert analyte_mass to analyte_volume if necessary
    if pd.isnull(analyte_mass):
        analyte_mass = (
            analyte_volume
            * density.seawater_1atm_MP81(
                temperature=temperature[0],
                salinity=salinity,
            )
            * 1e-3
        )
    return Converted(titrant_mass, dd.measurement, temperature, analyte_mass)


def get_totals_k_constants(converted, salinity, **kwargs):
    """Get total salt contents and equilibrium constants through a titration.

    Parameters
    ----------
    converted : Converted
        A namedtuple generated by `convert_amount_units`.
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
    totals, totals_pyco2 = interface.get_totals(salinity, **kwargs_totals)
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


def _get_solve_func(titrant, measurement_type):
    if titrant.upper() == "H2SO4":
        solve_func = core.solve_emf_complete_H2SO4
    else:
        if measurement_type.lower() == "emf":
            solve_func = core.solve_emf_complete
        else:
            solve_func = core.solve_emf_pH_adjust
    return solve_func


def calibrate(
    alkalinity_certified,
    converted,
    totals,
    k_constants,
    emf0_guess=None,
    measurement_type="emf",
    pH_min=3,
    pH_max=4,
    titrant_molinity_guess=0.1,
    titrant="HCl",
):
    """Solve for `titrant_molinity` given `alkalinity_certified`.

    Parameters
    ----------
    alkalinity_certified : float
        The target total alkalinity in µmol/kg-solution.
    converted : Converted
        The output from `convert_amount_units`.
    totals : dict
        The total salt contents through the titration, including dilution
        by the titrant, assuming the titrant is not one of the totals (e.g.,
        it is not H2SO4).  Generated with `get_totals_k_constants`.
    k_constants : dict
        The equilibrium constants through the titration, including dilution
        by the titrant (affects pH scale conversions only) and temperature
        variations.  Generated with `get_totals_k_constants`.
    emf0_guess : float, optional
        A first-guess value for EMF0 in mV.
    measurement_type : str, optional
        The type of measurement in the data file, "emf" (default) or "pH".
    pH_min : float, optional
        Minimum pH to use from the titration data, by default 3.
    pH_max : float, optional
        Maximum pH to use from the titration data, by default 4.
    titrant_molinity_guess : float, optional
        First guess for the molinity of titrant in mol/kg-solution, by default
        0.1.
    titrant : str, optional
        What the titrant was, "HCl" (default) or "H2SO4".

    Returns
    -------
    float
        The least-squares best fitting `titrant_molinity`.
    """
    cv = converted
    titrant_molinity = core.calibrate(
        alkalinity_certified,
        cv.titrant_mass,
        cv.measurement,
        cv.temperature,
        cv.analyte_mass,
        totals,
        k_constants,
        pH_min=pH_min,
        pH_max=pH_max,
        solve_func=_get_solve_func(titrant, measurement_type),
        titrant_molinity_guess=titrant_molinity_guess,
    )["x"][0]
    return titrant_molinity


def solve(
    titrant_molinity,
    converted,
    totals,
    k_constants,
    emf0_guess=None,
    measurement_type="emf",
    pH_min=3,
    pH_max=4,
    titrant="HCl",
):
    """Solve alkalinity, EMF0 and initial pH for a titration file given the
    `titrant_molinity`.

    Parameters
    ----------
    titrant_molinity : float
        The titrant molinity in mol/kg-solution.
    converted : Converted
        The output from `convert_amount_units`.
    totals : dict
        The total salt contents through the titration, including dilution
        by the titrant, assuming the titrant is not one of the totals (e.g.,
        it is not H2SO4).  Generated with `get_totals_k_constants`.
    k_constants : dict
        The equilibrium constants through the titration, including dilution
        by the titrant (affects pH scale conversions only) and temperature
        variations.  Generated with `get_totals_k_constants`.
    emf0_guess : float, optional
        A first-guess value for EMF0 in mV.
    measurement_type : str, optional
        The type of measurement in the data file, "emf" (default) or "pH".
    pH_min : float, optional
        Minimum pH to use from the titration data, by default 3.
    pH_max : float, optional
        Maximum pH to use from the titration data, by default 4.
    titrant : str, optional
        What the titrant was, "HCl" (default) or "H2SO4".

    Returns
    -------
    SolveResult - a named tuple containing the fields
        alkalinity : float
            Total alkalinity in µmol/kg-solution.
        emf0 : float
            The electrode EMF0 in mV.
        pH_initial : float
            pH at the first measurement point on the free scale.
        temperature_initial : float
            Temperature at the first measurement point in °C.
        analyte_mass : float
            Mass of the analyte in kg.
        opt_result : dict
            The alkalinity-EMF0 solver output (see docstring for the
            appropriate solver in `core`).
    """
    cv = converted
    solve_func = _get_solve_func(titrant, measurement_type)
    opt_result = solve_func(
        titrant_molinity,
        cv.titrant_mass,
        cv.measurement,
        cv.temperature,
        cv.analyte_mass,
        totals,
        k_constants,
        emf0_guess=emf0_guess,
        pH_min=pH_min,
        pH_max=pH_max,
    )
    alkalinity, emf0 = opt_result["x"]
    # Calculate initial pH
    pH_initial = convert.pH_to_emf(cv.measurement[0], 0, cv.temperature[0])
    pH_initial = convert.emf_to_pH(pH_initial, emf0, cv.temperature[0])
    return SolveResult(
        alkalinity * 1e6,
        emf0,
        pH_initial,
        cv.temperature[0],
        cv.analyte_mass,
        opt_result,
    )


def _get_kwarg_keys(func):
    params = inspect.signature(func).parameters
    return {p for p in params if params[p].default != inspect.Parameter.empty}


# These are the kwargs that can be passed to `get_totals_k_constants`
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
# These are the kwargs that can be passed to other functions
keys_cau = _get_kwarg_keys(convert_amount_units)
keys_calibrate = _get_kwarg_keys(calibrate)
keys_solve = _get_kwarg_keys(solve)
