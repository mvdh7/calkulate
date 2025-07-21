# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019--2025  Matthew P. Humphreys  (GNU GPLv3)
"""Work with datasets containing multiple titrations."""

from warnings import warn

import numpy as np
import pandas as pd
import PyCO2SYS as pyco2

from . import files
from .core import SolveEmfResult, SolvePhGranResult, SolvePhResult
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
            cal = files.calibrate(
                row.file_name,
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
        .apply(get_group_calibration, include_groups=False)
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
    ds = solve(ds, verbose=verbose, **kwargs)
    return ds


def add_solve_results(solved, sr):
    if isinstance(sr, SolveEmfResult):
        solved["alkalinity_npts"] = sr.used.sum()
        solved["alkalinity_std"] = sr.alkalinity_all[sr.used].std()
        solved["alkalinity"] = sr.alkalinity
        solved["emf0"] = sr.emf0
        solved["gran_alkalinity"] = sr.ggr.alkalinity * 1e6
        solved["gran_emf0"] = sr.ggr.emf0
        solved["pH_init"] = sr.pH[0]
        solved["temperature_init"] = sr.temperature[0]
    elif isinstance(sr, SolvePhResult):
        solved["alkalinity_npts"] = sr.used.sum()
        solved["alkalinity_std"] = sr.alkalinity_all[sr.used].std()
        solved["alkalinity"] = sr.alkalinity
        solved["pH_init"] = sr.pH[0]
        solved["temperature_init"] = sr.temperature[0]
    elif isinstance(sr, SolvePhGranResult):
        solved["alkalinity_npts"] = sr.used.sum()
        solved["alkalinity"] = sr.alkalinity
        solved["pH_init"] = sr.pH[0]
        solved["temperature_init"] = sr.temperature[0]
    return solved


def solve_row(row, verbose=False, **kwargs):
    """Solve alkalinity, EMF0 and initial pH for one titration in a dataset."""
    # Define blank output
    solved = pd.Series(
        {
            "alkalinity_npts": 0,
            "alkalinity_std": np.nan,
            "alkalinity": np.nan,
            "analyte_mass": row.analyte_mass,
            "emf0": np.nan,
            "gran_alkalinity": np.nan,
            "gran_emf0": np.nan,
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
            sr = files.solve(
                row.file_name,
                row.titrant_molinity,
                row.salinity,
                **kwargs_solve,
            )
            solved = add_solve_results(solved, sr)
        except Exception as e:
            print(f'Error solving "{row.file_name}":')
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
    # Check for bad kwargs, but don't break on them
    kwargs_ignored = []
    for k in _backcompat(kwargs.copy(), []):
        if k not in files.keys_solve | {"pH_range", "read_dat_kwargs"}:
            kwargs_ignored.append(k)
    if len(kwargs_ignored) > 0:
        warn(
            "kwargs not recognised, being ignored: "
            + ("{} " * len(kwargs_ignored)).format(*kwargs_ignored)
        )
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
