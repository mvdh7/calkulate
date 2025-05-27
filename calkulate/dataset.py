# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019--2025  Matthew P. Humphreys  (GNU GPLv3)
"""Work with datasets containing multiple titrations."""

import os

import numpy as np
import pandas as pd
import PyCO2SYS as pyco2

# from . import default, titration
from .read.titrations import keys_read_dat, read_dat
from .titration import (
    calibrate as t_calibrate,
    convert_amount_units,
    get_totals_k_constants,
    keys_calibrate,
    keys_cau,
    keys_solve,
    keys_totals_ks,
    solve as t_solve,
)


def get_total_salts(ds):
    """Estimate total salt contents from salinity using PyCO2SYS without
    overwriting existing values.  Operates in-place.
    """
    assert "salinity" in ds, 'Dataset must contain a "salinity" column.'
    # Use opt_total_borate = 1 where it's not provided
    if "opt_total_borate" in ds:
        ds["opt_total_borate"] = ds.opt_total_borate.where(
            ds.opt_total_borate.notnull(), 1
        )
    else:
        ds["opt_total_borate"] = 1
    results = pyco2.sys(
        salinity=ds.salinity.values,
        opt_total_borate=ds.opt_total_borate.values,
    )
    salts = ["total_sulfate", "total_borate", "total_fluoride"]
    for salt in salts:
        if salt not in ds:
            ds[salt] = np.nan
        ds[salt] = ds[salt].where(ds[salt].notnull(), other=results[salt])
    return ds


def _backcompat(row, kwargs):
    """Deal with old-style kwargs for backwards compatibility."""
    prevdicts = ["read_dat_kwargs"]
    for prevdict in prevdicts:
        if prevdict in kwargs:
            for k, v in kwargs[prevdict].items():
                kwargs[k] = v
    if "pH_range" in kwargs:
        kwargs["pH_min"], kwargs["pH_max"] = kwargs["pH_range"]
    old2new = {"read_dat_method": "file_type"}
    for old, new in old2new.items():
        if old in kwargs:
            kwargs[new] = kwargs[old]
        if "read_dat_method" in row:
            kwargs[new] = row[old]
    return row, kwargs


def _get_kwargs_for(keys, kwargs, row):
    # Start by getting any kwargs that are in the keys
    kwargs_for = {k: v for k, v in kwargs.items() if k in keys}
    # Overwrite any that have values in the row
    for k in keys:
        if k in row and pd.notnull(row[k]):
            kwargs_for[k] = row[k]
    return kwargs_for


def calibrate_row(row, verbose=False, **kwargs):
    """Calibrate `titrant_molinity` for a single row of (i.e., a single
    titration in) a dataset.  Also returns analyte mass.
    """
    row, kwargs = _backcompat(row, kwargs)
    # Initialise output
    calibrated = pd.Series(
        {
            "titrant_molinity_here": np.nan,
            "analyte_mass": row.analyte_mass,
        }
    )
    if pd.notnull(row.alkalinity_certified) and row.file_good:
        if verbose:
            print(f"Calkulate: calibrating {row.file_name}...")
        # STEP 1: IMPORT TITRATION DATA
        # Append file_name to file_path if present
        if "file_path" in row:
            file_name = os.path.join(row.file_path, row.file_name)
        else:
            file_name = row.file_name
        # Get kwargs for read_dat
        kwargs_read_dat = _get_kwargs_for(keys_read_dat, kwargs, row)
        # Import the file
        try:
            dat_data = read_dat(file_name, **kwargs_read_dat)
        except FileNotFoundError:
            print(f'File not found: "{file_name}"')
            return calibrated
        # STEP 2: CONVERT AMOUNT UNITS
        # Get kwargs for convert_amount_units
        kwargs_cau = _get_kwargs_for(keys_cau, kwargs, row)
        # If there is analyte_mass, use it, otherwise use analyte_volume
        if "analyte_volume" in kwargs_cau and "analyte_mass" in kwargs_cau:
            kwargs_cau["analyte_volume"] = None
        # Prepare titration file, totals and k_constants
        try:
            converted = convert_amount_units(
                dat_data, row.salinity, **kwargs_cau
            )
            calibrated["analyte_mass"] = converted.analyte_mass
        except Exception as e:
            print(f'Error converting units for "{file_name}":')
            print(f"{e}")
            return calibrated
        # STEP 3: GET TOTAL SALTS AND EQUILIBRIUM CONSTANTS
        kwargs_totals_ks = _get_kwargs_for(keys_totals_ks, kwargs, row)
        try:
            totals, k_constants = get_totals_k_constants(
                converted,
                row.salinity,
                **kwargs_totals_ks,
            )
        except Exception as e:
            print(f'Error getting salts and constants for "{file_name}":')
            print(f"{e}")
            return calibrated
        # STEP 4: DO THE CALIBRATION
        kwargs_calibrate = _get_kwargs_for(keys_calibrate, kwargs, row)
        try:
            calibrated["titrant_molinity_here"] = t_calibrate(
                row.alkalinity_certified,
                converted,
                totals,
                k_constants,
                **kwargs_calibrate,
            )
        except Exception as e:
            print(f'Error calibrating "{file_name}":')
            print(f"{e}")
            return calibrated
    return calibrated


def get_group_calibration(ds_group):
    """Get mean titrant molinity and statistics for each analysis_batch group."""
    titrant_molinities = ds_group.titrant_molinity_here[
        ds_group.reference_good & ~np.isnan(ds_group.titrant_molinity_here)
    ]
    return pd.Series(
        {
            "titrant_molinity": titrant_molinities.mean(),
            "titrant_molinity__std": titrant_molinities.std(),
            "titrant_molinity__count": np.size(titrant_molinities),
        }
    )


def get_batches(ds):
    """Get mean titrant molinity and statistics for all analysis_batch groups."""
    batches = (
        ds[["analysis_batch", "titrant_molinity_here", "reference_good"]]
        .groupby(by="analysis_batch")
        .apply(get_group_calibration)
    )
    return batches


def calibrate(ds, verbose=False, **kwargs):
    """Calibrate `titrant_molinity` for all titrations with an
    `alkalinity_certified` value and assign means based on `analysis_batch`.

    Parameters
    ----------
    ds : pandas.DataFrame
        A table containing metadata for each titration (not used if running as
        a method).
    verbose : bool, optional
        Whether to print progress, by default `calk.default.verbose`.

    Returns
    -------
    pandas.DataFrame
        The titration metadataset with additional columns found by the solver.
    """
    print("Calkulate: calibrating titrant_molinity...")
    # Get analyte_mass from analyte_volume if required
    if "analyte_mass" not in ds:
        assert "analyte_volume" in ds, (
            "ds must contain either 'analyte_mass' or 'analyte_volume'!"
        )
        ds["analyte_mass"] = np.nan
    # Check essential columns are present
    for must_have in ["alkalinity_certified", "analyte_mass", "salinity"]:
        assert must_have in ds, f'ds must contain a "{must_have}" column!'
    if "titrant_amount_unit" in ds:
        ds["titrant_amount_unit"] = np.where(
            pd.isnull(ds.titrant_amount_unit),
            "ml",
            ds.titrant_amount_unit,
        )
    # Calibrate titrant_molinity_here for each row with an alkalinity_certified
    if "file_good" not in ds:
        ds["file_good"] = True
    calibrated_rows = ds.apply(
        calibrate_row, axis=1, verbose=verbose, **kwargs
    )
    for k, v in calibrated_rows.items():
        ds[k] = v
    # Get titrant_molinity averaged by analysis_batch
    if "analysis_batch" not in ds:
        ds["analysis_batch"] = 0
    if "reference_good" not in ds:
        ds["reference_good"] = ~np.isnan(ds.titrant_molinity_here)
    batches = get_batches(ds)
    ds["titrant_molinity"] = batches.loc[
        ds.analysis_batch, "titrant_molinity"
    ].to_numpy()
    print("Calkulate: calibration complete!")
    return ds


def solve_row(row, verbose=False, **kwargs):
    """Solve alkalinity, EMF0 and initial pH for one titration in a dataset."""
    row, kwargs = _backcompat(row, kwargs)
    # Define blank output
    solved = pd.Series(
        {
            "alkalinity_guess": np.nan,
            "emf0_guess": np.nan,
            "alkalinity": np.nan,
            "alkalinity_std": np.nan,
            "alkalinity_npts": 0,
            "emf0": np.nan,
            "pH_initial": np.nan,
            "temperature_initial": np.nan,
            "analyte_mass": row.analyte_mass,
        }
    )
    if pd.notnull(row.titrant_molinity) and row.file_good:
        if verbose:
            print(f"Calkulate: solving {row.file_name}...")
        # STEP 1: IMPORT TITRATION DATA
        # Append file_name to file_path if present
        if "file_path" in row:
            file_name = os.path.join(row.file_path, row.file_name)
        else:
            file_name = row.file_name
        # Get kwargs for read_dat
        kwargs_read_dat = _get_kwargs_for(keys_read_dat, kwargs, row)
        # Import the file
        try:
            dat_data = read_dat(file_name, **kwargs_read_dat)
        except FileNotFoundError:
            print(f'File not found: "{file_name}"')
            return solved
        except Exception as e:
            print(f'Error importing "{file_name}":')
            print(f"{e}")
            return solved
        # STEP 2: CONVERT AMOUNT UNITS
        # Get kwargs for convert_amount_units
        kwargs_cau = _get_kwargs_for(keys_cau, kwargs, row)
        # If there is analyte_mass, use it, otherwise use analyte_volume
        if "analyte_volume" in kwargs_cau and "analyte_mass" in kwargs_cau:
            kwargs_cau["analyte_volume"] = None
        # Prepare titration file, totals and k_constants
        try:
            converted = convert_amount_units(
                dat_data, row.salinity, **kwargs_cau
            )
            solved["analyte_mass"] = converted.analyte_mass
        except Exception as e:
            print(f'Error converting units for "{file_name}":')
            print(f"{e}")
            return solved
        # STEP 3: GET TOTAL SALTS AND EQUILIBRIUM CONSTANTS
        kwargs_totals_ks = _get_kwargs_for(keys_totals_ks, kwargs, row)
        try:
            totals, k_constants = get_totals_k_constants(
                converted,
                row.salinity,
                **kwargs_totals_ks,
            )
        except Exception as e:
            print(f'Error getting salts and constants for "{file_name}":')
            print(f"{e}")
            return solved
        # STEP 4: SOLVE THE TITRATION
        kwargs_solve = _get_kwargs_for(keys_solve, kwargs, row)
        try:
            sr = t_solve(
                row.titrant_molinity,
                converted,
                totals,
                k_constants,
                **kwargs_solve,
            )
            solved["alkalinity_guess"] = sr.opt_result["alkalinity_guess"]
            solved["emf0_guess"] = sr.opt_result["emf0_guess"]
            solved["alkalinity"] = sr.alkalinity
            solved["alkalinity_std"] = np.std(sr.opt_result["fun"]) * 1e6
            solved["alkalinity_npts"] = sum(sr.opt_result["data_used"])
            solved["emf0"] = sr.emf0
            solved["pH_initial"] = sr.pH_initial
            solved["temperature_initial"] = sr.temperature_initial
        except Exception as e:
            print(f'Error calibrating "{file_name}":')
            print(f"{e}")
            return solved
    return solved


def solve(ds, verbose=False, **kwargs):
    """Solve alkalinity, EMF0 and initial pH for all titrations with a
    `titrant_molinity` value in a `Dataset`.

    Parameters
    ----------
    ds : pandas.DataFrame
        A table containing metadata for each titration.
    verbose : `bool`, optional
        Whether to print progress, by default False.

    Returns
    -------
    pandas.DataFrame
        The titration metadataset with additional columns found by the solver.
    """
    print("Calkulate: solving alkalinity...")
    # Get analyte_mass from analyte_volume if required
    if "analyte_mass" not in ds:
        assert "analyte_volume" in ds, (
            "ds must contain either 'analyte_mass' or 'analyte_volume'!"
        )
        ds["analyte_mass"] = np.nan
    if "file_good" not in ds:
        ds["file_good"] = True
    if "titrant_amount_unit" in ds:
        ds["titrant_amount_unit"] = np.where(
            pd.isnull(ds.titrant_amount_unit),
            "ml",
            ds.titrant_amount_unit,
        )
    solved_rows = ds.apply(solve_row, axis=1, verbose=verbose, **kwargs)
    for k, v in solved_rows.items():
        ds[k] = v
    if "alkalinity_certified" in ds:
        ds["alkalinity_offset"] = ds.alkalinity - ds.alkalinity_certified
    print("Calkulate: solving complete!")
    return ds


def calkulate(ds, verbose=False, **kwargs):
    """Calibrate and then solve all titrations in a `Dataset`.

    Parameters
    ----------
    ds : pd.DataFrame
        A table containing metadata for each titration (not used if running as a
        method).
    verbose : `bool`, optional
        Whether to print progress, by default `calk.default.verbose`.

    Returns
    -------
    pd.DataFrame
        The titration metadataset with additional columns found by the solver.
    """
    get_total_salts(ds)
    calibrate(ds, verbose=verbose, **kwargs)
    solve(ds, verbose=verbose, **kwargs)
    return ds


def to_Titration(ds, index, read_dat_kwargs={}):
    """Create a `calk.Titration` object from one row of a `Dataset`.

    Parameters
    ----------
    ds : `pandas.DataFrame` or `calk.Dataset`
        The `Dataset` to make the Titration from (not used if running as a method).
    index
        The row index in the `Dataset` to use.
    read_dat_kwargs : `dict`, optional
        Any kwargs that need to be passed through in order to correctly read the
        titration data file (e.g., `encoding`).

    Returns
    -------
    `calk.Titration`
        A `calk.Titration` for the specified row of the `Dataset`.
    """
    dsr = ds.loc[index]
    prepare_kwargs = {"read_dat_kwargs": read_dat_kwargs}
    for k, v in prepare_defaults.items():
        if k in dsr:
            if not pd.isnull(dsr[k]):
                prepare_kwargs[k] = dsr[k]
            else:
                prepare_kwargs[k] = v
        else:
            prepare_kwargs[k] = v
    analyte_mass = prepare_kwargs.pop("analyte_mass")
    analyte_volume = prepare_kwargs.pop("analyte_volume")
    tt = titration.Titration(
        file_name=dsr.file_name,
        file_path=dsr.file_path if "file_path" in dsr else "",
        salinity=dsr.salinity,
        analyte_mass=analyte_mass,
        analyte_volume=analyte_volume,
        file_prepare_kwargs=prepare_kwargs,
    )
    if "alkalinity_certified" in dsr:
        if not pd.isnull(dsr.alkalinity_certified):
            tt.alkalinity_certified = dsr.alkalinity_certified
    if "titrant_molinity" in dsr:
        if not pd.isnull(dsr.titrant_molinity):
            tt.set_titrant_molinity(dsr.titrant_molinity)
            tt.solve()
    return tt
