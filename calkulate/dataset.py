# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019--2025  Matthew P. Humphreys  (GNU GPLv3)
"""Work with datasets containing multiple titrations."""

from collections import namedtuple

import numpy as np
import pandas as pd
import PyCO2SYS as pyco2

from . import default, titration


def get_total_salts(ds):
    """Estimate total salt contents from salinity using PyCO2SYS without
    overwriting existing values.  Operates in-place.
    """
    assert "salinity" in ds, "ds must contain a 'salinity' column!"
    if "opt_total_borate" in ds:
        opt_total_borate = ds.opt_total_borate.where(
            ds.opt_total_borate.notnull(), default.opt_total_borate
        ).to_numpy()
    else:
        opt_total_borate = 1
    results = pyco2.sys(
        salinity=ds.salinity.to_numpy(),
        opt_total_borate=opt_total_borate,
    )
    salts = ["total_sulfate", "total_borate", "total_fluoride"]
    for salt in salts:
        if salt not in ds:
            ds[salt] = np.nan
        ds[salt] = ds[salt].where(ds[salt].notnull(), other=results[salt])
    return ds


RowArgs = namedtuple(
    "RowArgs",
    ("file_name",),
)

# These are the kwargs that can be passed to `titration.get_totals_k_constants`
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


def calibrate_row(row, kwargs_dat_data=None, verbose=False):
    """Calibrate `titrant_molinity` for a single row of (i.e., a single
    titration in) a dataset."""
    if pd.notnull(row.alkalinity_certified) & row.file_good:
        if verbose:
            print(f"Calkulate: calibrating {row.file_name}...")
        # Append file_name to file_path if present
        if "file_path" in row:
            file_name = row.file_path + row.file_name
        else:
            file_name = row.file_name
        # If there is analyte_mass, use it, otherwise use analyte_volume
        analyte_amount = {}
        if "analyte_mass" in row and pd.notnull(row.analyte_mass):
            analyte_amount["analyte_mass"] = row.analyte_mass
        else:
            analyte_amount["analyte_volume"] = row.analyte_volume
        # Get row kwargs for get_dat_data
        if kwargs_dat_data is None:
            kwargs_dat_data = {}
        if "temperature_override" in row and pd.notnull(
            row.temperature_override
        ):
            kwargs_dat_data["temperature_override"] = row.temperature_override
        # Prepare titration file, totals and k_constants
        try:
            prepped = titration.prepare(
                file_name,
                row.salinity,
                **analyte_amount,
                kwargs_dat_data=kwargs_dat_data,
            )
            analyte_mass = prepped.analyte_mass
        except FileNotFoundError:
            print(f'Calkulate - file not found: "{row.file_name}"')
            return pd.Series(
                {
                    "titrant_molinity_here": np.nan,
                    "analyte_mass": row.analyte_mass,
                }
            )
        kwargs_totals_ks = {}
        for k in keys_totals_ks:
            if k in row and pd.notnull(row[k]):
                kwargs_totals_ks[k] = row[k]
        totals, k_constants = titration.get_totals_k_constants(
            prepped.titrant_mass,
            prepped.temperature,
            prepped.analyte_mass,
            row.salinity,
            **kwargs_totals_ks,
        )
        # Do the calibration
        kwargs_calibrate = {}
        for k in [
            "emf0_guess",
            "measurement_type",
            "pH_min",
            "pH_max",
            "titrant_molinity_guess",
            "titrant",
        ]:
            if k in row and pd.notnull(row[k]):
                kwargs_calibrate[k] = row[k]
        try:
            titrant_molinity_here = titration.calibrate(
                row.alkalinity_certified,
                prepped,
                totals,
                k_constants,
                **kwargs_calibrate,
            )
        except Exception as e:
            print(f"Calkulate ERROR calibrating '{row.file_name}': {e}")
            titrant_molinity_here = np.nan
            analyte_mass = row.analyte_mass
    else:
        # If alkalinity_certified not provided for this row
        titrant_molinity_here = np.nan
        analyte_mass = row.analyte_mass
    return pd.Series(
        {
            "titrant_molinity_here": titrant_molinity_here,
            "analyte_mass": analyte_mass,
        }
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


def calibrate(ds, kwargs_dat_data=None, verbose=False):
    """Calibrate `titrant_molinity` for all titrations with an
    `alkalinity_certified` value and assign means based on `analysis_batch`.

    Parameters
    ----------
    ds : pandas.DataFrame or calk.Dataset
        A table containing metadata for each titration (not used if running as
        a method).
    kwargs_dat_data : dict, optional
        kwargs to pass to `get_dat_data`, by default `None`.
    verbose : bool, optional
        Whether to print progress, by default `calk.default.verbose`.

    Returns
    -------
    pandas.DataFrame or calk.Dataset
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
            default.titrant_amount_unit,
            ds.titrant_amount_unit,
        )
    # Calibrate titrant_molinity_here for each row with an alkalinity_certified
    if "file_good" not in ds:
        ds["file_good"] = True
    calibrated_rows = ds.apply(
        calibrate_row,
        axis=1,
        kwargs_dat_data=kwargs_dat_data,
        verbose=verbose,
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


def solve_row(row, kwargs_dat_data=None, verbose=False):
    """Solve alkalinity, EMF0 and initial pH for one titration in a dataset."""
    if verbose:
        print("Calkulate: solving {}...".format(row.file_name))
    solved = False
    if pd.notnull(row.titrant_molinity) & row.file_good:
        # Append file_name to file_path if present
        if "file_path" in row:
            file_name = row.file_path + row.file_name
        else:
            file_name = row.file_name
        # If there is analyte_mass, use it, otherwise use analyte_volume
        analyte_amount = {}
        if "analyte_mass" in row and pd.notnull(row.analyte_mass):
            analyte_amount["analyte_mass"] = row.analyte_mass
        else:
            analyte_amount["analyte_volume"] = row.analyte_volume
        # Get row kwargs for get_dat_data
        if kwargs_dat_data is None:
            kwargs_dat_data = {}
        if "temperature_override" in row and pd.notnull(
            row.temperature_override
        ):
            kwargs_dat_data["temperature_override"] = row.temperature_override
        # Prepare titration file, totals and k_constants
        try:
            prepped = titration.prepare(
                file_name,
                row.salinity,
                **analyte_amount,
                kwargs_dat_data=kwargs_dat_data,
            )
        except FileNotFoundError:
            print(f'Calkulate - file not found: "{row.file_name}"')
            return pd.Series(
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
        kwargs_totals_ks = {}
        for k in keys_totals_ks:
            if k in row and pd.notnull(row[k]):
                kwargs_totals_ks[k] = row[k]
        totals, k_constants = titration.get_totals_k_constants(
            prepped.titrant_mass,
            prepped.temperature,
            prepped.analyte_mass,
            row.salinity,
            **kwargs_totals_ks,
        )
        # TODO up to here was similar to `calibrate_row`; make a separate
        # subfunction for both?
        kwargs_solve = {}
        for k in [
            "emf0_guess",
            "measurement_type",
            "pH_min",
            "pH_max",
            "titrant",
        ]:
            if k in row and pd.notnull(row[k]):
                kwargs_solve[k] = row[k]
        try:
            solved = titration.solve(
                row.titrant_molinity,
                prepped,
                totals,
                k_constants,
                **kwargs_solve,
            )
        except Exception as e:
            print(f'Calkulate ERROR solving "{row.file_name}":  {e}')
    if solved:
        return pd.Series(
            {
                "alkalinity_guess": solved.opt_result["alkalinity_guess"],
                "emf0_guess": solved.opt_result["emf0_guess"],
                "alkalinity": solved.alkalinity,
                "alkalinity_std": np.std(solved.opt_result["fun"]) * 1e6,
                "alkalinity_npts": sum(solved.opt_result["data_used"]),
                "emf0": solved.emf0,
                "pH_initial": solved.pH_initial,
                "temperature_initial": solved.temperature_initial,
                "analyte_mass": solved.analyte_mass,
            }
        )
    else:
        return pd.Series(
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


def solve(ds, kwargs_dat_data=None, verbose=False):
    """Solve alkalinity, EMF0 and initial pH for all titrations with a
    `titrant_molinity` value in a `Dataset`.

    Parameters
    ----------
    ds : `pandas.DataFrame` or `calk.Dataset`
        A table containing metadata for each titration (not used if running as
        a method).
    kwargs_dat_data : dict, optional
        kwargs to pass to `get_dat_data`, by default `None`.
    verbose : `bool`, optional
        Whether to print progress, by default `calk.default.verbose`.

    Returns
    -------
    pandas.DataFrame or calk.Dataset
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
            default.titrant_amount_unit,
            ds.titrant_amount_unit,
        )
    solved_rows = ds.apply(
        solve_row,
        axis=1,
        kwargs_dat_data=kwargs_dat_data,
        verbose=verbose,
    )
    for k, v in solved_rows.items():
        ds[k] = v
    if "alkalinity_certified" in ds:
        ds["alkalinity_offset"] = ds.alkalinity - ds.alkalinity_certified
    print("Calkulate: solving complete!")
    return ds


def calkulate(ds, kwargs_dat_data=None, verbose=default.verbose):
    """Calibrate and then solve all titrations in a `Dataset`.

    Parameters
    ----------
    ds : `pandas.DataFrame` or `calk.Dataset`
        A table containing metadata for each titration (not used if running as a
        method).
    kwargs_dat_data : dict, optional
        kwargs to pass to `get_dat_data`, by default `None`.
    verbose : `bool`, optional
        Whether to print progress, by default `calk.default.verbose`.

    Returns
    -------
    `pandas.DataFrame` or `calk.Dataset`
        The titration metadataset with additional columns found by the solver.
    """
    get_total_salts(ds)
    calibrate(ds, kwargs_dat_data=kwargs_dat_data, verbose=verbose)
    solve(ds, kwargs_dat_data=kwargs_dat_data, verbose=verbose)
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


class Dataset(pd.DataFrame):
    """
    `calk.Dataset`
    ================
    A `calk.Dataset` is pandas `DataFrame` with several `calk.dataset` functions
    available as methods.

    Alkalinity solving methods
    --------------------------
    `calibrate`
        Find the best-fit `titrant_molinity` for each sample that has a value for
        `alkalinity_certified`.
    `solve`
        Solve every sample for `alkalinity` when `titrant_molinity` is known.
    `calkulate`
        Run the `calibrate` and `solve` steps sequentially.

    Data visualisation methods
    --------------------------
    `plot_titrant_molinity`
        Plot the `titrant_molinity` values determined with `calibrate` through time.
    `plot_alkalinity_offset`
        Plot the difference between solved `alkalinity` and `alkalinity_certified`
        through time.

    Conversion methods
    ------------------
    `to_Titration`
        Return a `calk.Titration` object for one row in the `Dataset`.
    `to_pandas`
        Return a copy of the `Dataset` as a standard pandas `DataFrame`.
    """

    get_batches = get_batches
    get_total_salts = get_total_salts
    calibrate = calibrate
    solve = solve
    calkulate = calkulate
    to_Titration = to_Titration

    # from .plot import (
    #     alkalinity_offset as plot_alkalinity_offset,
    #     titrant_molinity as plot_titrant_molinity,
    # )

    def to_pandas(self):
        """Return a copy of the `Dataset` as a standard pandas `DataFrame`."""
        return pd.DataFrame(self)
