# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019--2025  Matthew P. Humphreys  (GNU GPLv3)
"""Work with titration data in a file.

This is the middle level of the three layers of processing functions.
It includes convenience functions for importing data files and reformatting
them to use with the `core` functions, one file at a time.

The processing steps are
  1.  `read.read_dat`: import titration data file.
  2.  `convert.amount_units`: convert amount units to kg (titrant and analyte).
  3.  `core.totals_ks`: calculate total salt contents and equilibrium
      constants.
  4a. `calibrate`: find best fitting titrant molinity given a certified
      alkalinity value.
  4b. `solve`: find best fitting alkalinity and EMF0 given titrant molinity.
"""

from collections import namedtuple

from . import convert, core
from .meta import _get_kwarg_keys


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
        The target total alkalinity in µmol/kg-sol.
    converted : Converted
        The output from `convert.amount_units`.
    totals : dict
        The total salt contents through the titration, including dilution
        by the titrant, assuming the titrant is not one of the totals (e.g.,
        it is not H2SO4).  Generated with `totals_ks`.
    k_constants : dict
        The equilibrium constants through the titration, including dilution
        by the titrant (affects pH scale conversions only) and temperature
        variations.  Generated with `totals_ks`.
    emf0_guess : float, optional
        A first-guess value for EMF0 in mV.
    measurement_type : str, optional
        The type of measurement in the data file, "emf" (default) or "pH".
    pH_min : float, optional
        Minimum pH to use from the titration data, by default 3.
    pH_max : float, optional
        Maximum pH to use from the titration data, by default 4.
    titrant_molinity_guess : float, optional
        First guess for the molinity of titrant in mol/kg-sol, by default 0.1.
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
        The titrant molinity in mol/kg-sol.
    converted : Converted
        The output from `convert.amount_units`.
    totals : dict
        The total salt contents through the titration, including dilution
        by the titrant, assuming the titrant is not one of the totals (e.g.,
        it is not H2SO4).  Generated with `totals_ks`.
    k_constants : dict
        The equilibrium constants through the titration, including dilution
        by the titrant (affects pH scale conversions only) and temperature
        variations.  Generated with `totals_ks`.
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
            Total alkalinity in µmol/kg-sol.
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


keys_calibrate = _get_kwarg_keys(calibrate)
keys_solve = _get_kwarg_keys(solve)
