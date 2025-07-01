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
  4a. `core.calibrate_*`: find best fitting titrant molinity given a certified
      alkalinity value.
  4b. `core.solve_*`: find best fitting alkalinity and EMF0 given titrant
      molinity.
"""

import os
from warnings import warn

from . import core
from .convert import amount_units, keys_cau, pH_to_emf
from .core import (
    add_titrant_totals,
    keys_calibrate_emf,
    keys_calibrate_pH,
    keys_calibrate_pH_gran,
    keys_solve_emf,
    keys_solve_pH,
    keys_solve_pH_gran,
    keys_titrant_totals,
    keys_totals_ks,
    totals_ks,
)
from .meta import _get_kwargs_for
from .read.titrations import keys_read_dat, read_dat


keys_calibrate = (
    keys_read_dat
    | keys_cau
    | keys_totals_ks
    | keys_calibrate_emf
    | keys_calibrate_pH
    | keys_calibrate_pH_gran
    | {"solve_mode", "file_path"}
)
keys_solve = (
    keys_read_dat
    | keys_cau
    | keys_totals_ks
    | keys_titrant_totals
    | keys_solve_emf
    | keys_solve_pH
    | keys_solve_pH_gran
    | {"solve_mode", "file_path"}
)


def calibrate(
    file_name,
    alkalinity_certified,
    salinity,
    solve_mode="emf",
    **kwargs,
):
    """Solve for `titrant_molinity` given `alkalinity_certified`.

    Parameters
    ----------
    file_name : str
        The name (and path to) the titration data file.
    alkalinity_certified : float
        The target total alkalinity in µmol/kg-sol.
    salinity : float
        Practical salinity of the analyte.
    solve_mode : str, optional, case-insensitive
        How to solve for alkalinity, one of
            "emf" (default) - measurements are EMF in mV
            "pH_adjust" - measurements are pH but their EMF0 can be adjusted
            "pH" - measurements are pH and cannot be adjusted
            "pH_gran" - measurements are pH, use Gran-plot solver
    kwargs
        Any keyword arguments that need passing to lower-level functions
        (`read_dat`, `amount_units`, `totals_ks` and `calibrate_*`).

    Returns
    -------
    opt_result : scipy.optimize.OptimizeResult
        Output from `scipy.optimize.least_squares`, where the solved value is
        `titrant_molinity = opt_result["x"][0]`.
    """
    # Check for bad kwargs, but don't break on them
    kwargs_ignored = []
    for k in kwargs:
        if k not in keys_calibrate:
            kwargs_ignored.append(k)
    if len(kwargs_ignored) > 0:
        warn(
            "kwargs not recognised, being ignored: "
            + ("{} " * len(kwargs_ignored)).format(*kwargs_ignored)
        )
    # Import the titration data file
    if "file_path" in kwargs:
        file_name = os.path.join(kwargs["file_path"], file_name)
    kwargs_read_dat = _get_kwargs_for(keys_read_dat, kwargs)
    dd = read_dat(file_name, **kwargs_read_dat)
    # Convert amount units
    kwargs_cau = _get_kwargs_for(keys_cau, kwargs)
    cv = amount_units(dd, salinity, **kwargs_cau)
    # Get total salts and equilibrium constants
    kwargs_totals_ks = _get_kwargs_for(keys_totals_ks, kwargs)
    totals, k_constants = totals_ks(cv, **kwargs_totals_ks)
    # Calibrate!
    if solve_mode.lower() == "emf":
        # Titration data are EMFs
        kwargs_calibrate_emf = _get_kwargs_for(keys_calibrate_emf, kwargs)
        cal = core.calibrate_emf(
            alkalinity_certified,
            cv.titrant_mass,
            cv.measurement,
            cv.temperature,
            cv.analyte_mass,
            totals,
            k_constants,
            **kwargs_calibrate_emf,
        )
    elif solve_mode.lower() == "ph_adjust":
        # Titration data are pHs but we want to allow the EMF0 to be adjusted
        kwargs_calibrate_emf = _get_kwargs_for(keys_calibrate_emf, kwargs)
        emf0_init = 0
        kwargs_calibrate_emf["emf0_init"] = emf0_init
        emf = pH_to_emf(cv.measurement, emf0_init, cv.temperature)
        cal = core.calibrate_emf(
            alkalinity_certified,
            cv.titrant_mass,
            emf,
            cv.temperature,
            cv.analyte_mass,
            totals,
            k_constants,
            **kwargs_calibrate_emf,
        )
    elif solve_mode.lower() == "ph":
        # Titration data are pHs and cannot be adjusted
        kwargs_calibrate_pH = _get_kwargs_for(keys_calibrate_pH, kwargs)
        cal = core.calibrate_pH(
            alkalinity_certified,
            cv.titrant_mass,
            cv.measurement,
            cv.temperature,
            cv.analyte_mass,
            totals,
            k_constants,
            **kwargs_calibrate_pH,
        )
    elif solve_mode.lower() == "ph_gran":
        # Titration data are pHs and cannot be adjusted
        kwargs_calibrate_pH_gran = _get_kwargs_for(
            keys_calibrate_pH_gran, kwargs
        )
        cal = core.calibrate_pH_gran(
            alkalinity_certified,
            cv.titrant_mass,
            cv.measurement,
            cv.temperature,
            cv.analyte_mass,
            totals,
            k_constants,
            **kwargs_calibrate_pH_gran,
        )
    else:
        raise Exception("`solve_mode` not valid")
    return cal


def solve(
    file_name,
    titrant_molinity,
    salinity,
    solve_mode="emf",
    **kwargs,
):
    """Solve for `alkalinity` etc. given `titrant_molinity`.

    Parameters
    ----------
    file_name : str
        The name (and path to) the titration data file.
    titrant_molinity : float
        The titrant molinity in mol/kg-sol.
    salinity : float
        Practical salinity of the analyte.
    solve_mode : str, optional, case-insensitive
        How to solve for alkalinity, one of
            "emf" (default) - measurements are EMF in mV
            "pH_adjust" - measurements are pH but their EMF0 can be adjusted
            "pH" - measurements are pH and cannot be adjusted
    kwargs
        Any keyword arguments that need passing to lower-level functions
        (`read_dat`, `amount_units`, `totals_ks`, `add_titrant_totals` and
        `solve_*`).

    Returns
    -------
    Depends on solve_mode:

    """
    # Import the titration data file
    if "file_path" in kwargs:
        file_name = os.path.join(kwargs["file_path"], file_name)
    kwargs_read_dat = _get_kwargs_for(keys_read_dat, kwargs)
    dd = read_dat(file_name, **kwargs_read_dat)
    # Convert amount units
    kwargs_cau = _get_kwargs_for(keys_cau, kwargs)
    cv = amount_units(dd, salinity, **kwargs_cau)
    # Get total salts and equilibrium constants
    kwargs_totals_ks = _get_kwargs_for(keys_totals_ks, kwargs)
    totals, k_constants = totals_ks(cv, **kwargs_totals_ks)
    kwargs_titrant_totals = _get_kwargs_for(keys_titrant_totals, kwargs)
    totals = add_titrant_totals(
        totals,
        cv.titrant_mass,
        cv.analyte_mass,
        titrant_molinity,
        titrant_molinity_prev=0,
        **kwargs_titrant_totals,
    )
    # Calibrate!
    if solve_mode.lower() == "emf":
        # Titration data are EMFs
        kwargs_solve_emf = _get_kwargs_for(keys_solve_emf, kwargs)
        sr = core.solve_emf(
            titrant_molinity,
            cv.titrant_mass,
            cv.measurement,
            cv.temperature,
            cv.analyte_mass,
            totals,
            k_constants,
            **kwargs_solve_emf,
        )
    elif solve_mode.lower() == "ph_adjust":
        # Titration data are pHs but we want to allow the EMF0 to be adjusted
        kwargs_solve_emf = _get_kwargs_for(keys_solve_emf, kwargs)
        emf0_init = 0
        kwargs_solve_emf["emf0_init"] = emf0_init
        emf = pH_to_emf(cv.measurement, emf0_init, cv.temperature)
        sr = core.solve_emf(
            titrant_molinity,
            cv.titrant_mass,
            emf,
            cv.temperature,
            cv.analyte_mass,
            totals,
            k_constants,
            **kwargs_solve_emf,
        )
    elif solve_mode.lower() == "ph":
        # Titration data are pHs and cannot be adjusted
        kwargs_solve_pH = _get_kwargs_for(keys_solve_pH, kwargs)
        sr = core.solve_pH(
            titrant_molinity,
            cv.titrant_mass,
            cv.measurement,
            cv.temperature,
            cv.analyte_mass,
            totals,
            k_constants,
            **kwargs_solve_pH,
        )
    elif solve_mode.lower() == "ph_gran":
        # Titration data are pHs and cannot be adjusted
        kwargs_solve_pH_gran = _get_kwargs_for(keys_solve_pH_gran, kwargs)
        sr = core.solve_pH_gran(
            titrant_molinity,
            cv.titrant_mass,
            cv.measurement,
            cv.temperature,
            cv.analyte_mass,
            totals,
            k_constants,
            **kwargs_solve_pH_gran,
        )
    else:
        raise Exception("`solve_mode` not valid")
    return sr
