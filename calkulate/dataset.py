# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019--2021  Matthew P. Humphreys  (GNU GPLv3)
"""Work with datasets containing multiple titrations."""

import copy
import numpy as np, pandas as pd
from . import default, titration


prepare_defaults = dict(
    analyte_mass=None,  # kg
    analyte_volume=None,  # ml
    dic=0,
    total_alpha=0,
    total_beta=0,
    total_ammonia=0,
    total_phosphate=0,
    total_silicate=0,
    total_sulfide=0,
    total_borate=None,
    total_fluoride=None,
    total_sulfate=None,
    k_alpha=None,
    k_ammonia=None,
    k_beta=None,
    k_bisulfate=None,
    k_borate=None,
    k_carbonic_1=None,
    k_carbonic_2=None,
    k_fluoride=None,
    k_phosphoric_1=None,
    k_phosphoric_2=None,
    k_phosphoric_3=None,
    k_silicate=None,
    k_sulfide=None,
    k_water=None,
    molinity_HCl=default.molinity_HCl,
    molinity_NaCl=default.molinity_NaCl,
    temperature_override=None,
    titrant_amount_unit=default.titrant_amount_unit,
    opt_k_bisulfate=default.opt_k_bisulfate,
    opt_k_carbonic=default.opt_k_carbonic,
    opt_k_fluoride=default.opt_k_fluoride,
    opt_total_borate=default.opt_total_borate,
    read_dat_method=default.read_dat_method,
)


def get_prepare_kwargs(ds_row):
    """Get prepare_kwargs for calibrate_row and solve_row functions."""
    # Get kwargs for titration.prepare
    prepare_kwargs = {}
    for k in prepare_defaults:
        if k in ds_row:
            if ~pd.isnull(ds_row[k]):
                prepare_kwargs[k] = ds_row[k]
    # Add analyte_mass or analyte_volume as needed
    if ~pd.isnull(ds_row.analyte_mass):
        prepare_kwargs["analyte_mass"] = ds_row.analyte_mass
    else:
        prepare_kwargs["analyte_volume"] = ds_row.analyte_volume
    return prepare_kwargs


def calibrate_row(
    ds_row,
    pH_range=default.pH_range,
    least_squares_kwargs=default.least_squares_kwargs,
    read_dat_kwargs={},
    verbose=default.verbose,
):
    """Calibrate titrant_molinity for all titrations with an alkalinity_certified
    value and assign means based on analysis_batch.
    """
    if ~np.isnan(ds_row.alkalinity_certified) & ds_row.file_good:
        if verbose:
            print("Calkulate: calibrating {}...".format(ds_row.file_name))
        prepare_kwargs = get_prepare_kwargs(ds_row)
        prepare_kwargs["read_dat_kwargs"] = read_dat_kwargs
        # Calibrate titrant molinity
        if "file_path" in ds_row:
            file_name = ds_row.file_path + ds_row.file_name
        else:
            file_name = ds_row.file_name
        titrant_molinity_guess = default.titrant_molinity_guess
        if "titrant_molinity_guess" in ds_row:
            if ~pd.isnull(ds_row.titrant_molinity_guess):
                titrant_molinity_guess = ds_row.titrant_molinity_guess
        try:
            titrant_molinity_here, analyte_mass = titration.calibrate(
                file_name,
                ds_row.salinity,
                ds_row.alkalinity_certified,
                titrant_molinity_guess=titrant_molinity_guess,
                pH_range=pH_range,
                least_squares_kwargs=least_squares_kwargs,
                **prepare_kwargs,
            )
        except FileNotFoundError:
            print("Calkulate: file not found: '{}'".format(ds_row.file_name))
            titrant_molinity_here = np.nan
            analyte_mass = ds_row.analyte_mass
        except:
            print("Calkulate: ERROR calibrating '{}'!".format(ds_row.file_name))
            titrant_molinity_here = np.nan
            analyte_mass = ds_row.analyte_mass
    else:
        # If alkalinity_certified not provided for this ds_row
        titrant_molinity_here = np.nan
        analyte_mass = ds_row.analyte_mass
    return pd.Series(
        {"titrant_molinity_here": titrant_molinity_here, "analyte_mass": analyte_mass,}
    )


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


def calibrate(
    ds,
    least_squares_kwargs=default.least_squares_kwargs,
    pH_range=default.pH_range,
    read_dat_kwargs={},
    inplace=True,
    verbose=default.verbose,
):
    """Calibrate titrant_molinity for all titrations with an alkalinity_certified
    value and assign means based on analysis_batch.
    """
    print("Calkulate: calibrating titrant_molinity...")
    # Get analyte_mass from analyte_volume if required
    if "analyte_mass" not in ds:
        assert (
            "analyte_volume" in ds
        ), "ds must contain either 'analyte_mass' or 'analyte_volume'!"
        ds["analyte_mass"] = np.nan
    # Check essential columns are present
    for must_have in ["alkalinity_certified", "analyte_mass", "salinity"]:
        assert must_have in ds, "ds must contain a '{}' column!".format(must_have)
    if not inplace:
        ds = copy.deepcopy(ds)
    if "titrant_amount_unit" in ds:
        ds["titrant_amount_unit"] = np.where(
            pd.isnull(ds.titrant_amount_unit),
            default.titrant_amount_unit,
            ds.titrant_amount_unit,
        )
    # Calibrate titrant_molinity_here for each row with an alkalinity_certified
    if "file_good" not in ds:
        ds["file_good"] = True
    calibrated_rows = ds.apply(
        calibrate_row,
        axis=1,
        pH_range=pH_range,
        least_squares_kwargs=least_squares_kwargs,
        read_dat_kwargs=read_dat_kwargs,
        verbose=verbose,
    )
    for k, v in calibrated_rows.iteritems():
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


def solve_row(
    ds_row,
    pH_range=default.pH_range,
    least_squares_kwargs=default.least_squares_kwargs,
    read_dat_kwargs={},
    verbose=default.verbose,
):
    """Solve alkalinity, EMF0 and initial pH for one titration in a dataset."""
    if verbose:
        print("Calkulate: solving {}...".format(ds_row.file_name))
    if ~np.isnan(ds_row.titrant_molinity) & ds_row.file_good:
        prepare_kwargs = get_prepare_kwargs(ds_row)
        prepare_kwargs["read_dat_kwargs"] = read_dat_kwargs
        # Solve for alkalinity etc.
        if "file_path" in ds_row:
            file_name = ds_row.file_path + ds_row.file_name
        else:
            file_name = ds_row.file_name
        try:
            (
                alkalinity,
                emf0,
                pH_initial,
                temperature_initial,
                analyte_mass,
            ) = titration.solve(
                file_name,
                ds_row.salinity,
                ds_row.titrant_molinity,
                pH_range=pH_range,
                least_squares_kwargs=least_squares_kwargs,
                **prepare_kwargs,
            )
        except FileNotFoundError:
            print("Calkulate: file not found: '{}'".format(ds_row.file_name))
            alkalinity = emf0 = pH_initial = temperature_initial = np.nan
            analyte_mass = ds_row.analyte_mass
        except:
            print("Calkulate: ERROR solving '{}'!".format(ds_row.file_name))
            alkalinity = emf0 = pH_initial = temperature_initial = np.nan
            analyte_mass = ds_row.analyte_mass
    else:
        # If alkalinity_certified not provided for this ds_row
        alkalinity = emf0 = pH_initial = temperature_initial = np.nan
        analyte_mass = ds_row.analyte_mass
    return pd.Series(
        {
            "alkalinity": alkalinity,
            "emf0": emf0,
            "pH_initial": pH_initial,
            "temperature_initial": temperature_initial,
            "analyte_mass": analyte_mass,
        }
    )


def solve(
    ds,
    least_squares_kwargs=default.least_squares_kwargs,
    pH_range=default.pH_range,
    read_dat_kwargs={},
    inplace=True,
    verbose=default.verbose,
):
    """Solve alkalinity, EMF0 and initial pH for all titrations with a
    titrant_molinity value in a dataset.
    """
    print("Calkulate: solving alkalinity...")
    if not inplace:
        ds = copy.deepcopy(ds)
    if "file_good" not in ds:
        ds["file_good"] = True
    if "titrant_amount_unit" in ds:
        ds["titrant_amount_unit"] = np.where(
            pd.isnull(ds.titrant_amount_unit),
            default.titrant_amount_unit,
            ds.titrant_amount_unit,
        )
    solved_rows = ds.apply(
        solve_row,
        axis=1,
        least_squares_kwargs=least_squares_kwargs,
        pH_range=pH_range,
        read_dat_kwargs=read_dat_kwargs,
        verbose=verbose,
    )
    for k, v in solved_rows.iteritems():
        ds[k] = v
    print("Calkulate: solving complete!")
    return ds


def calkulate(
    ds,
    least_squares_kwargs=default.least_squares_kwargs,
    pH_range=default.pH_range,
    read_dat_kwargs={},
    inplace=True,
    verbose=default.verbose,
):
    """Calibrate and solve all titrations in a dataset."""
    ds = calibrate(
        ds,
        pH_range=pH_range,
        least_squares_kwargs=least_squares_kwargs,
        read_dat_kwargs=read_dat_kwargs,
        inplace=inplace,
        verbose=verbose,
    )
    ds = solve(
        ds,
        pH_range=pH_range,
        least_squares_kwargs=least_squares_kwargs,
        read_dat_kwargs=read_dat_kwargs,
        inplace=inplace,
        verbose=verbose,
    )
    return ds


class Dataset(pd.DataFrame):
    """pandas DataFrame with dataset functions available as methods."""

    get_batches = get_batches
    calibrate = calibrate
    solve = solve
    calkulate = calkulate

    from .plot import (
        titrant_molinity as plot_titrant_molinity,
        alkalinity_offset as plot_alkalinity_offset,
    )
