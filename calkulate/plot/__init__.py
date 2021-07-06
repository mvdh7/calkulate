# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019--2021  Matthew P. Humphreys  (GNU GPLv3)
"""Visualise the results."""

import numpy as np
from matplotlib import pyplot as plt
from .. import meta, dataset
from . import titration, misc


def titrant_molinity(
    data,
    xvar=None,
    show_bad=True,
    show_batches=True,
    figure_fname=None,
):
    """Plot the individually calibrated titrant_molinity values and batch averages."""
    if xvar is None:
        xdata = np.arange(len(data.index))
        xlabel = "Analysis number"
    else:
        xdata = data[xvar]
        xlabel = xvar
    if "reference_good" in data:
        G = data.reference_good.values
    else:
        G = np.full(np.size(xdata), True)
    # Draw the figure
    fig, ax = plt.subplots(dpi=300)
    ax.scatter(
        xdata[G],
        data.titrant_molinity_here[G],
        alpha=0.6,
        c="xkcd:navy",
        s=20,
        label="Good data",
    )
    if show_bad:
        ax.scatter(
            xdata[~G],
            data.titrant_molinity_here[~G],
            alpha=0.8,
            c="xkcd:strawberry",
            s=20,
            marker="x",
            label="Bad data",
        )
        ax.legend()
    if show_batches:
        batches = dataset.get_batches(data)
        for batch in batches.index:
            B = (data.analysis_batch == batch).to_numpy()
            ax.plot(xdata[B], data.titrant_molinity[B], c="xkcd:navy")
    ax.grid(alpha=0.3)
    misc.add_credit(ax)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"Titrant molinity / mol$\cdot$kg$^{-1}$")
    plt.tight_layout()
    if figure_fname is not None:
        plt.savefig(figure_fname)
    return fig, ax


def alkalinity_offset(
    data,
    xvar=None,
    show_bad=True,
    show_batches=True,
    figure_fname=None,
):
    """Plot the offset between measured and certified values for reference materials
    after solving with batch-averaged titrant_molinity.
    """
    if xvar is None:
        xdata = np.arange(len(data.index))
        xlabel = "Analysis number"
    else:
        xdata = data[xvar]
        xlabel = xvar
    if "reference_good" in data:
        G = data.reference_good.values
    else:
        G = np.full(np.size(xdata), True)
    if show_batches:
        batches = dataset.get_batches(data)
        clr = np.full(np.size(xdata), -1)
        for b, batch in enumerate(batches.index):
            B = (data.analysis_batch == batch).to_numpy()
            clr[B] = b
        clr = clr[G]
    else:
        clr = "xkcd:navy"
    # Draw the figure
    fig, ax = plt.subplots(dpi=300)
    ax.scatter(
        xdata[G],
        data.alkalinity_offset[G],
        alpha=0.6,
        c=clr,
        cmap="turbo",
        s=20,
        label="Good data",
    )
    if show_bad:
        ax.scatter(
            xdata[~G],
            data.alkalinity_offset[~G],
            alpha=0.8,
            c="xkcd:strawberry",
            s=20,
            marker="x",
            label="Bad data",
        )
        ax.legend()
    ax.grid(alpha=0.3)
    ax.axhline(0, c="k", lw=0.8)
    misc.add_credit(ax)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$\Delta$ Alkalinity (measured $-$ certified) / mol$\cdot$kg$^{-1}$")
    plt.tight_layout()
    if figure_fname is not None:
        plt.savefig(figure_fname)
    return fig, ax
