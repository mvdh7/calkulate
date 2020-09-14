import pandas as pd, numpy as np
from .. import options


def prepare(dataset, read_dat_kwargs=None):
    """Do everything except for calibrating and solving."""
    if read_dat_kwargs is None:
        read_dat_kwargs = {}
    dataset.get_titrations(**read_dat_kwargs)
    dataset.get_analyte_mass()
    dataset.get_titrant_mass()
    dataset.get_totals()
    dataset.get_k_constants()
    return dataset


def calibrate(dataset, index=None, pH_range=(3, 4), verbose=options.verbose):
    """Calibrate titrant molinity for all samples with known alkalinity."""
    if index is None:
        dataset.calibrate_all(pH_range=pH_range, verbose=verbose)
    dataset.set_batch_mean_molinity()
    return dataset


def solve(dataset, index=None, pH_range=(3, 4), verbose=options.verbose):
    """Solve all samples with known titrant molinity for total alkalinity."""
    if index is None:
        dataset.solve_all(pH_range=pH_range, verbose=verbose)
    return dataset


def set_batch_mean_molinity(dataset):
    """Calculate mean titrant molinity for each analysis batch and apply to all samples
    in that batch.
    """
    if "analysis_batch" not in dataset:
        dataset["analysis_batch"] = 0
    use_titrations = ~np.isnan(dataset.titrant_molinity_here)
    if "reference_good" in dataset:
        use_titrations &= dataset["reference_good"]
    batches = (
        dataset[["analysis_batch", "titrant_molinity_here"]][use_titrations]
        .groupby("analysis_batch")
        .agg(["mean", "std", "count"])
    )
    dataset["titrant_molinity"] = batches.loc[
        dataset.analysis_batch.values, "titrant_molinity_here"
    ]["mean"].values
    return dataset


def calkulate(dataset, read_dat_kwargs=None, pH_range=(3, 4), verbose=options.verbose):
    """Do absolutely everything in one go."""
    dataset.prepare(read_dat_kwargs=read_dat_kwargs)
    dataset.calibrate_all(pH_range=pH_range, verbose=verbose)
    dataset.set_batch_mean_molinity()
    dataset.solve_all(pH_range=pH_range, verbose=verbose)
    return dataset


class Dataset(pd.DataFrame):
    from .get import (
        get_measurement_type,
        get_titrations,
        get_analyte_temperature,
        get_analyte_mass,
        get_titrant_density,
        get_titrant_mass,
        get_analyte_totals,
        get_titration_totals,
        get_totals,
        get_k_constants,
    )
    from .quantify import calibrate_all, solve_all

    prepare = prepare
    calibrate = calibrate
    set_batch_mean_molinity = set_batch_mean_molinity
    solve = solve
    calkulate = calkulate
