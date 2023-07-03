# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019--2023  Matthew P. Humphreys  (GNU GPLv3)
"""Work with datasets containing multiple titrations."""

import copy
import numpy as np, pandas as pd
import PyCO2SYS as pyco2
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
    molinity_H2SO4=default.molinity_H2SO4,
    temperature_override=None,
    titrant_amount_unit=default.titrant_amount_unit,
    titrant_density=None,
    opt_k_bisulfate=default.opt_k_bisulfate,
    opt_k_carbonic=default.opt_k_carbonic,
    opt_k_fluoride=default.opt_k_fluoride,
    # opt_pH_scale=default.opt_pH_scale,  # doesn't work yet
    opt_total_borate=default.opt_total_borate,
    read_dat_method=default.read_dat_method,
)


def get_total_salts(ds, inplace=True):
    """Estimate total salt contents from salinity using PyCO2SYS, without overwriting
    existing values.
    """
    assert "salinity" in ds, "ds must contain a 'salinity' column!"
    if "opt_total_borate" in ds:
        opt_total_borate = ds.opt_total_borate.where(
            ~pd.isnull(ds.opt_total_borate), default.opt_total_borate
        ).to_numpy()
    else:
        opt_total_borate = default.opt_total_borate
    results = pyco2.sys(
        8,
        2000,
        3,
        2,
        salinity=ds.salinity.to_numpy(),
        opt_total_borate=opt_total_borate,
    )
    if not inplace:
        ds = copy.deepcopy(ds)
    salts = ["total_sulfate", "total_borate", "total_fluoride"]
    for salt in salts:
        if salt not in ds:
            ds[salt] = np.nan
        ds[salt].where(~pd.isnull(ds[salt]), other=results[salt], inplace=True)
    return ds


def get_prepare_kwargs(ds_row):
    """Get prepare_kwargs for calibrate_row and solve_row functions."""
    # Get kwargs for titration.prepare
    prepare_kwargs = {}
    for k in prepare_defaults:
        if k in ds_row:
            if not pd.isnull(ds_row[k]):
                prepare_kwargs[k] = ds_row[k]
    # Add analyte_mass or analyte_volume as needed
    if not pd.isnull(ds_row.analyte_mass):
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
            if not pd.isnull(ds_row.titrant_molinity_guess):
                titrant_molinity_guess = ds_row.titrant_molinity_guess
        emf0_guess = None
        if "emf0_guess" in ds_row:
            if not pd.isnull(ds_row.emf0_guess):
                emf0_guess = ds_row.emf0_guess
        # Deal with H2SO4 titrant special case
        titrant = default.titrant
        analyte_total_sulfate = None
        if "titrant" in ds_row:
            if not pd.isnull(ds_row.titrant):
                titrant = ds_row.titrant
                if titrant == "H2SO4":
                    assert "total_sulfate" in ds_row
                    assert not pd.isnull(ds_row.total_sulfate)
                    analyte_total_sulfate = ds_row.total_sulfate
        # Calibrate!
        try:
            titrant_molinity_here, analyte_mass = titration.calibrate(
                file_name,
                ds_row.salinity,
                ds_row.alkalinity_certified,
                analyte_total_sulfate=analyte_total_sulfate,
                titrant=titrant,
                titrant_molinity_guess=titrant_molinity_guess,
                pH_range=pH_range,
                least_squares_kwargs=least_squares_kwargs,
                emf0_guess=emf0_guess,
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


def calibrate(
    ds,
    least_squares_kwargs=default.least_squares_kwargs,
    pH_range=default.pH_range,
    read_dat_kwargs={},
    inplace=True,
    verbose=default.verbose,
):
    """Calibrate ``titrant_molinity`` for all titrations with an
    ``alkalinity_certified`` value and assign means based on ``analysis_batch``.

    Parameters
    ----------
    ds : ``pandas.DataFrame`` or ``calk.Dataset``
        A table containing metadata for each titration (not used if running as a
        method).
    least_squares_kwargs : ``dict``, optional
        kwargs to pass to scipy.optimize.least_squares, by default
        `calk.default.least_squares_kwargs`.
    pH_range : ``list``, optional
        The minimum and maximum pH values to use to solve, by default
        ``calk.default.pH_range``.
    read_dat_kwargs : ``dict``, optional
        kwargs to pass to ``read_dat``, by default ``{}``.
    inplace : ``bool``, optional
        Whether to perform the operation in-place on ``ds``, by default ``True``
    verbose : ``bool``, optional
        Whether to print progress, by default ``calk.default.verbose``.

    Returns
    -------
    ``pandas.DataFrame`` or ``calk.Dataset``
        The titration metadataset with additional columns found by the solver.
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
    solved = False
    if ~np.isnan(ds_row.titrant_molinity) & ds_row.file_good:
        prepare_kwargs = get_prepare_kwargs(ds_row)
        prepare_kwargs["read_dat_kwargs"] = read_dat_kwargs
        # Solve for alkalinity etc.
        if "file_path" in ds_row:
            file_name = ds_row.file_path + ds_row.file_name
        else:
            file_name = ds_row.file_name
        titrant = default.titrant
        analyte_total_sulfate = None
        if "titrant" in ds_row:
            if not pd.isnull(ds_row.titrant):
                titrant = ds_row.titrant
                if titrant == "H2SO4":
                    assert "total_sulfate" in ds_row
                    assert not pd.isnull(ds_row.total_sulfate)
                    analyte_total_sulfate = ds_row.total_sulfate
        emf0_guess = None
        if "emf0_guess" in ds_row:
            if not pd.isnull(ds_row.emf0_guess):
                emf0_guess = ds_row.emf0_guess
        try:
            (
                alkalinity,
                emf0,
                pH_initial,
                temperature_initial,
                analyte_mass,
                opt_result,
            ) = titration.solve(
                file_name,
                ds_row.salinity,
                ds_row.titrant_molinity,
                analyte_total_sulfate=analyte_total_sulfate,
                titrant=titrant,
                pH_range=pH_range,
                least_squares_kwargs=least_squares_kwargs,
                emf0_guess=emf0_guess,
                **prepare_kwargs,
            )
            solved = True
        except FileNotFoundError:
            print("Calkulate: file not found: '{}'".format(ds_row.file_name))
        except:
            print("Calkulate: ERROR solving '{}'!".format(ds_row.file_name))
    if solved:
        return pd.Series(
            {
                "alkalinity_gran": opt_result["alkalinity_gran"],
                "emf0_gran": opt_result["emf0_gran"],
                "alkalinity": alkalinity,
                "alkalinity_std": np.std(opt_result["fun"]) * 1e6,
                "alkalinity_npts": sum(opt_result["data_used"]),
                "emf0": emf0,
                "pH_initial": pH_initial,
                "temperature_initial": temperature_initial,
                "analyte_mass": analyte_mass,
            }
        )
    else:
        return pd.Series(
            {
                "alkalinity_gran_estimate": np.nan,
                "emf0_gran_estimate": np.nan,
                "alkalinity": np.nan,
                "alkalinity_std": np.nan,
                "alkalinity_npts": 0,
                "emf0": np.nan,
                "pH_initial": np.nan,
                "temperature_initial": np.nan,
                "analyte_mass": ds_row.analyte_mass,
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
    ``titrant_molinity`` value in a ``Dataset``.

    Parameters
    ----------
    ds : ``pandas.DataFrame`` or ``calk.Dataset``
        A table containing metadata for each titration (not used if running as a
        method).
    least_squares_kwargs : ``dict``, optional
        kwargs to pass to scipy.optimize.least_squares, by default
        `calk.default.least_squares_kwargs`.
    pH_range : ``list``, optional
        The minimum and maximum pH values to use to solve, by default
        ``calk.default.pH_range``.
    read_dat_kwargs : ``dict``, optional
        kwargs to pass to ``read_dat``, by default ``{}``.
    inplace : ``bool``, optional
        Whether to perform the operation in-place on ``ds``, by default ``True``
    verbose : ``bool``, optional
        Whether to print progress, by default ``calk.default.verbose``.

    Returns
    -------
    ``pandas.DataFrame`` or ``calk.Dataset``
        The titration metadataset with additional columns found by the solver.
    """
    print("Calkulate: solving alkalinity...")
    if not inplace:
        ds = copy.deepcopy(ds)
    # Get analyte_mass from analyte_volume if required
    if "analyte_mass" not in ds:
        assert (
            "analyte_volume" in ds
        ), "ds must contain either 'analyte_mass' or 'analyte_volume'!"
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
        least_squares_kwargs=least_squares_kwargs,
        pH_range=pH_range,
        read_dat_kwargs=read_dat_kwargs,
        verbose=verbose,
    )
    for k, v in solved_rows.items():
        ds[k] = v
    if "alkalinity_certified" in ds:
        ds["alkalinity_offset"] = ds.alkalinity - ds.alkalinity_certified
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
    """Calibrate and then solve all titrations in a ``Dataset``.

    Parameters
    ----------
    ds : ``pandas.DataFrame`` or ``calk.Dataset``
        A table containing metadata for each titration (not used if running as a
        method).
    least_squares_kwargs : ``dict``, optional
        kwargs to pass to scipy.optimize.least_squares, by default
        `calk.default.least_squares_kwargs`.
    pH_range : ``list``, optional
        The minimum and maximum pH values to use to solve, by default
        ``calk.default.pH_range``.
    read_dat_kwargs : ``dict``, optional
        kwargs to pass to ``read_dat``, by default ``{}``.
    inplace : ``bool``, optional
        Whether to perform the operation in-place on ``ds``, by default ``True``
    verbose : ``bool``, optional
        Whether to print progress, by default ``calk.default.verbose``.

    Returns
    -------
    ``pandas.DataFrame`` or ``calk.Dataset``
        The titration metadataset with additional columns found by the solver.
    """
    ds = get_total_salts(ds, inplace=inplace)
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


def to_Titration(ds, index, read_dat_kwargs={}):
    """Create a ``calk.Titration`` object from one row of a ``Dataset``.

    Parameters
    ----------
    ds : ``pandas.DataFrame`` or ``calk.Dataset``
        The ``Dataset`` to make the Titration from (not used if running as a method).
    index
        The row index in the ``Dataset`` to use.
    read_dat_kwargs : ``dict``, optional
        Any kwargs that need to be passed through in order to correctly read the
        titration data file (e.g., ``encoding``).

    Returns
    -------
    ``calk.Titration``
        A ``calk.Titration`` for the specified row of the ``Dataset``.
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
    ``calk.Dataset``
    ================
    A ``calk.Dataset`` is pandas ``DataFrame`` with several ``calk.dataset`` functions
    available as methods.

    Alkalinity solving methods
    --------------------------
    ``calibrate``
        Find the best-fit ``titrant_molinity`` for each sample that has a value for
        ``alkalinity_certified``.
    ``solve``
        Solve every sample for ``alkalinity`` when ``titrant_molinity`` is known.
    ``calkulate``
        Run the ``calibrate`` and ``solve`` steps sequentially.

    Data visualisation methods
    --------------------------
    ``plot_titrant_molinity``
        Plot the ``titrant_molinity`` values determined with ``calibrate`` through time.
    ``plot_alkalinity_offset``
        Plot the difference between solved ``alkalinity`` and ``alkalinity_certified``
        through time.

    Conversion methods
    ------------------
    ``to_Titration``
        Return a ``calk.Titration`` object for one row in the ``Dataset``.
    ``to_pandas``
        Return a copy of the ``Dataset`` as a standard pandas ``DataFrame``.
    """

    get_batches = get_batches
    get_total_salts = get_total_salts
    calibrate = calibrate
    solve = solve
    calkulate = calkulate
    to_Titration = to_Titration

    from .plot import (
        titrant_molinity as plot_titrant_molinity,
        alkalinity_offset as plot_alkalinity_offset,
    )

    def to_pandas(self):
        """Return a copy of the ``Dataset`` as a standard pandas ``DataFrame``."""
        return pd.DataFrame(self)
