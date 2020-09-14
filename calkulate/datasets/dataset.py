import pandas as pd


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


def calibrate(dataset, index=None):
    if index is None:
        dataset.calibrate_all()
    return dataset


def solve(dataset, index=None):
    if index is None:
        dataset.solve_all()
    return dataset


def calkulate(dataset, read_dat_kwargs=None):
    """Do absolutely everything in one go."""
    dataset.prepare()
    dataset.calibrate_all()
    dataset.solve_all()
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
    solve = solve
    calkulate = calkulate
