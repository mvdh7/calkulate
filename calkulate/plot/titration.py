import numpy as np
from matplotlib import pyplot as plt
from .. import default
from . import misc


alpha_gran = 0.8
color_gran = "xkcd:azure"
marker_gran = "s"
alpha_final = 0.6
color_final = "k"


def emf(tt, ax=None, **scatter_kwargs):
    """EMF changes as acid is added throughout a titration."""
    ttt = tt.titration
    if ax is None:
        fig, ax = plt.subplots(dpi=default.dpi)
    ax.scatter(ttt.titrant_mass * 1e3, ttt.emf, **scatter_kwargs)
    ax.set_xlabel("Titrant mass / g")
    ax.set_ylabel("EMF / mV")
    misc.add_credit(ax)
    return ax


def pH(tt, ax=None, show_gran=True, **scatter_kwargs):
    """pH changes as acid is added throughout a titration."""
    ttt = tt.titration
    if ax is None:
        fig, ax = plt.subplots(dpi=default.dpi)
    if show_gran:
        ax.scatter(
            ttt.titrant_mass[~ttt.G_gran] * 1e3,
            ttt.pH_gran[~ttt.G_gran],
            label="Gran — unused",
            c="none",
            edgecolor=color_gran,
            alpha=alpha_gran,
            marker=marker_gran,
        )
        ax.scatter(
            ttt.titrant_mass[ttt.G_gran] * 1e3,
            ttt.pH_gran[ttt.G_gran],
            label="Gran — used",
            c=color_gran,
            edgecolor=color_gran,
            alpha=alpha_gran,
            marker=marker_gran,
        )
    ax.scatter(
        ttt.titrant_mass[~ttt.G_final] * 1e3,
        ttt.pH[~ttt.G_final],
        label="Final — unused",
        c="none",
        edgecolor=color_final,
        alpha=alpha_final,
    )
    ax.scatter(
        ttt.titrant_mass[ttt.G_final] * 1e3,
        ttt.pH[ttt.G_final],
        label="Final — used",
        c=color_final,
        edgecolor=color_final,
        alpha=alpha_final,
    )
    pH_range = (
        tt.calibrate_kwargs["pH_range"]
        if "pH_range" in tt.calibrate_kwargs
        else default.pH_range
    )
    xlim = ax.get_xlim()
    ax.fill_between(
        xlim,
        *pH_range,
        color=color_final,
        alpha=0.2,
        edgecolor="none",
        zorder=-1,
        label="Target pH range"
    )
    ax.set_xlim(xlim)
    ax.set_xlabel("Titrant mass / g")
    ax.set_ylabel("pH (Free scale)")
    ax.legend()
    misc.add_credit(ax)
    return ax


def gran_emf0(tt, ax=None):
    """Plot Gran estimate of EMF0."""
    ttt = tt.titration
    if ax is None:
        fig, ax = plt.subplots(dpi=default.dpi)
    ax.scatter(
        ttt.titrant_mass[~ttt.G_gran] * 1e3,
        ttt.emf0_gran[~ttt.G_gran],
        label="Unused",
        c="none",
        edgecolor=color_gran,
        alpha=alpha_gran,
        marker=marker_gran,
    )
    ax.scatter(
        ttt.titrant_mass[ttt.G_gran] * 1e3,
        ttt.emf0_gran[ttt.G_gran],
        label="Used",
        c=color_gran,
        edgecolor=color_gran,
        alpha=alpha_gran,
        marker=marker_gran,
    )
    ax.set_xlabel("Titrant mass / g")
    ax.set_ylabel("EMF° estimate / mV")
    ax.axhline(tt.emf0_gran, c=color_gran, label="Gran EMF°")
    ax.legend()
    ax.set_title("Gran EMF° = {:.2f} mV".format(tt.emf0_gran))
    misc.add_credit(ax)
    return ax


def gran_alkalinity(tt, ax=None):
    """Plot Gran estimate of alkalinity."""
    ttt = tt.titration
    if ax is None:
        fig, ax = plt.subplots(dpi=default.dpi)
    ax.scatter(
        ttt.titrant_mass[~ttt.G_gran] * 1e3,
        ttt.gran_estimates[~ttt.G_gran],
        color="none",
        edgecolor=color_gran,
        alpha=alpha_gran,
        marker=marker_gran,
        label="Unused",
    )
    ax.scatter(
        ttt.titrant_mass[ttt.G_gran] * 1e3,
        ttt.gran_estimates[ttt.G_gran],
        color=color_gran,
        edgecolor=color_gran,
        alpha=alpha_gran,
        marker=marker_gran,
        label="Used",
    )
    ax.axhline(0, c="k", lw=0.8, zorder=-1)
    x_line = np.array([-tt.gran_intercept / tt.gran_slope, ttt.titrant_mass.max()])
    y_line = x_line * tt.gran_slope + tt.gran_intercept
    ax.plot(x_line * 1e3, y_line, color=color_gran, label="Fit")
    ax.axvline(
        x_line[0] * 1e3, color=color_gran, lw=1, dashes=[3, 3], label="Intercept"
    )
    ax.set_xlabel("Titrant mass / g")
    ax.set_ylabel("Gran function")
    ax.set_title("Gran alkalinity = {:.1f} μmol/kg-sw".format(tt.alkalinity_gran))
    ax.legend()
    misc.add_credit(ax)
    return ax
