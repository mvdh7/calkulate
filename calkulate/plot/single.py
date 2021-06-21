# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019--2021  Matthew P. Humphreys  (GNU GPLv3)
"""Make plots for single titrations."""

import pandas as pd
from matplotlib import pyplot as plt
from .. import core, default, titration


def _emf(ax, titrant_mass, emf, scatter_kwargs={}):
    """EMF change as acid is added throughout a titration."""
    ax.scatter(titrant_mass * 1e3, emf, **scatter_kwargs)
    ax.set_xlabel("Titrant mass / g")
    ax.set_ylabel("EMF / mV")
    return ax


def emf(ds_row, ax=None, scatter_kwargs={}, read_dat_kwargs={}):
    """EMF change as acid is added throughout a titration."""

    def in_row(field, value_default):
        value = value_default
        if field in ds_row:
            value = ds_row[field] if not pd.isnull(ds_row[field]) else value_default
        else:
            value = value_default
        return value

    # Get defaults
    file_path = in_row("file_path", "")
    titrant = in_row("titrant", default.titrant)
    titrant_amount_unit = in_row("titrant_amount_unit", default.titrant_amount_unit)
    titrant_density = in_row("titrant_density", None)
    read_dat_method = in_row("read_dat_method", default.read_dat_method)
    # Read .dat file
    titrant_mass, emf = titration.get_dat_data(
        file_path + ds_row.file_name,
        titrant=titrant,
        titrant_amount_unit=titrant_amount_unit,
        titrant_density=titrant_density,
        read_dat_method=read_dat_method,
        read_dat_kwargs={},
    )[:2]
    # Draw plot
    if ax is None:
        fig, ax = plt.subplots(dpi=default.dpi)
    ax = _emf(ax, titrant_mass, emf, scatter_kwargs=scatter_kwargs)
    return ax


def _gran_alkalinity(
    ax,
    titrant_mass,
    emf,
    temperature,
    analyte_mass,
    titrant_molinity,
    titrant=default.titrant,
    scatter_kwargs={},
):
    """Gran function for first alkalinity guess."""
    mixture_mass = titrant_mass + analyte_mass
    gran_estimates = core.gran_estimator(mixture_mass, emf, temperature)
    G = core.get_gran_G(gran_estimates)
    alkalinity_guess = core.gran_guess_alkalinity(
        titrant_mass[G],
        gran_estimates[G],
        analyte_mass,
        titrant_molinity,
    )[0]
    if titrant == "H2SO4":
        alkalinity_guess *= 2
    ax.scatter(titrant_mass, gran_estimates[~G], label="Unused")
    ax.scatter(titrant_mass, gran_estimates[G], label="Used")
    gradient, intercept_y = linregress(titrant_mass[G], gran_estimates[G])[:2]
    intercept_x = -intercept_y / gradient
    plot_x = np.array([intercept_x, np.max(gran_estimates[G])])
    plot_y = plot_x * gradient + intercept_y
    ax.plot(plot_x, plot_y)
    ax.axhline(0, c="k", lw=0.8)
    ax.legend()
    return ax
