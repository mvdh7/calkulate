# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019--2025  Matthew P. Humphreys  (GNU GPLv3)
"""Work with datasets containing multiple titrations."""

import os

import numpy as np
import pandas as pd
import PyCO2SYS as pyco2

from . import files
from .core import SolveEmfResult, SolvePhResult
from .meta import _get_kwargs_for


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
    salts = ["total_sulfate", "total_borate", "total_fluoride"]
    do_pyco2 = not np.all([s in ds and ds[s].notnull().all() for s in salts])
    if do_pyco2:
        results = pyco2.sys(
            salinity=ds.salinity.values,
            opt_total_borate=ds.opt_total_borate.values,
        )
        for salt in salts:
            if salt not in ds:
                ds[salt] = np.nan
            ds[salt] = ds[salt].where(ds[salt].notnull(), other=results[salt])
    return ds


def _backcompat(kwargs, row):
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
    return kwargs


def calibrate_row(row, verbose=False, **kwargs):
    """Calibrate `titrant_molinity` for a single row of (i.e., a single
    titration in) a dataset.
    """
    # Initialise output
    titrant_molinity_here = np.nan
    if pd.notnull(row.alkalinity_certified) and row.file_good:
        if verbose:
            print(f"Calibrating {row.file_name}...")
        try:
            kwargs = _backcompat(kwargs, row)
            kwargs = _get_kwargs_for(files.keys_calibrate, kwargs, row)
            if "file_path" in row and pd.notnull(row.file_path):
                file_name = os.path.join(row.file_path, row.file_name)
            elif "file_path" in kwargs and pd.notnull(kwargs["file_path"]):
                file_name = os.path.join(kwargs["file_path"], row.file_name)
            else:
                file_name = row.file_name
            cal = files.calibrate(
                file_name,
                row.alkalinity_certified,
                row.salinity,
                **kwargs,
            )
            titrant_molinity_here = cal["x"][0]
        except Exception as e:
            print(f'Error calibrating "{row.file_name}":')
            print(f"{e}")
            return titrant_molinity_here
    return titrant_molinity_here


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


def prepare(ds):
    get_total_salts(ds)
    # Get analyte_mass from analyte_volume if required
    if "analyte_mass" not in ds:
        assert "analyte_volume" in ds, (
            "ds must contain either 'analyte_mass' or 'analyte_volume'!"
        )
        ds["analyte_mass"] = np.nan
    # Check essential columns are present
    assert "salinity" in ds, 'ds must contain a "salinity" column!'
    if "titrant_amount_unit" in ds:
        ds["titrant_amount_unit"] = np.where(
            pd.isnull(ds.titrant_amount_unit),
            "ml",
            ds.titrant_amount_unit,
        )
    if "file_good" not in ds:
        ds["file_good"] = True


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
    prepare(ds)
    assert "alkalinity_certified" in ds, (
        'ds must contain an "alkalinity_certified" column!'
    )
    # Calibrate titrant_molinity_here for each row with an alkalinity_certified
    ds["titrant_molinity_here"] = ds.apply(
        calibrate_row, axis=1, verbose=verbose, **kwargs
    )
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
    # Define blank output
    solved = pd.Series(
        {
            "alkalinity_gran": np.nan,
            "alkalinity_npts": 0,
            "alkalinity_std": np.nan,
            "alkalinity": np.nan,
            "analyte_mass": row.analyte_mass,
            "emf0_gran": np.nan,
            "emf0": np.nan,
            "pH_init": np.nan,
            "temperature_init": np.nan,
        }
    )
    if pd.notnull(row.titrant_molinity) and row.file_good:
        if verbose:
            print(f"Solving {row.file_name}...")
        try:
            kwargs = _backcompat(kwargs, row)
            kwargs_solve = _get_kwargs_for(files.keys_solve, kwargs, row)
            if "file_path" in row and pd.notnull(row.file_path):
                file_name = os.path.join(row.file_path, row.file_name)
            elif "file_path" in kwargs and pd.notnull(kwargs["file_path"]):
                file_name = os.path.join(kwargs["file_path"], row.file_name)
            else:
                file_name = row.file_name
            sr = files.solve(
                file_name,
                row.titrant_molinity,
                row.salinity,
                **kwargs_solve,
            )
            if isinstance(sr, SolveEmfResult):
                solved["alkalinity_gran"] = sr.ggr.alkalinity * 1e6
                solved["alkalinity_npts"] = sr.used.sum()
                solved["alkalinity_std"] = sr.alkalinity_all[sr.used].std()
                solved["alkalinity"] = sr.alkalinity
                solved["emf0_gran"] = sr.ggr.emf0
                solved["emf0"] = sr.emf0
                solved["pH_init"] = sr.pH[0]
                solved["temperature_init"] = sr.temperature[0]
            elif isinstance(sr, SolvePhResult):
                solved["alkalinity_npts"] = sr.used.sum()
                solved["alkalinity_std"] = sr.alkalinity_all[sr.used].std()
                solved["alkalinity"] = sr.alkalinity
                solved["pH_init"] = sr.pH[0]
                solved["temperature_init"] = sr.temperature[0]
        except Exception as e:
            print(f'Error solving "{file_name}":')
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
    prepare(ds)
    assert "titrant_molinity" in ds, (
        'ds must contain an "titrant_molinity" column!'
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
    calibrate(ds, verbose=verbose, **kwargs)
    solve(ds, verbose=verbose, **kwargs)
    return ds


# def to_Titration(ds, index, read_dat_kwargs={}):
#     """Create a `calk.Titration` object from one row of a `Dataset`.

#     Parameters
#     ----------
#     ds : `pandas.DataFrame` or `calk.Dataset`
#         The `Dataset` to make the Titration from (not used if running as a method).
#     index
#         The row index in the `Dataset` to use.
#     read_dat_kwargs : `dict`, optional
#         Any kwargs that need to be passed through in order to correctly read the
#         titration data file (e.g., `encoding`).

#     Returns
#     -------
#     `calk.Titration`
#         A `calk.Titration` for the specified row of the `Dataset`.
#     """
#     dsr = ds.loc[index]
#     prepare_kwargs = {"read_dat_kwargs": read_dat_kwargs}
#     for k, v in prepare_defaults.items():
#         if k in dsr:
#             if not pd.isnull(dsr[k]):
#                 prepare_kwargs[k] = dsr[k]
#             else:
#                 prepare_kwargs[k] = v
#         else:
#             prepare_kwargs[k] = v
#     analyte_mass = prepare_kwargs.pop("analyte_mass")
#     analyte_volume = prepare_kwargs.pop("analyte_volume")
#     tt = titration.Titration(
#         file_name=dsr.file_name,
#         file_path=dsr.file_path if "file_path" in dsr else "",
#         salinity=dsr.salinity,
#         analyte_mass=analyte_mass,
#         analyte_volume=analyte_volume,
#         file_prepare_kwargs=prepare_kwargs,
#     )
#     if "alkalinity_certified" in dsr:
#         if not pd.isnull(dsr.alkalinity_certified):
#             tt.alkalinity_certified = dsr.alkalinity_certified
#     if "titrant_molinity" in dsr:
#         if not pd.isnull(dsr.titrant_molinity):
#             tt.set_titrant_molinity(dsr.titrant_molinity)
#             tt.solve()
#     return tt
