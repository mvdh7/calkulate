import numpy as np
from scipy.interpolate import pchip_interpolate
from matplotlib import pyplot as plt
from .. import default, simulate
from . import misc


class Styler:
    def __init__(self, alpha=0.8, color="k", marker="o"):
        self.alpha = alpha
        self.color = color
        self.marker = marker

    def scatter(self, used=True):
        kwargs = {
            "alpha": self.alpha,
            "c": self.color if used else "none",
            "edgecolor": self.color,
            "marker": self.marker,
        }
        return kwargs

    def plot(self):
        kwargs = {
            "color": self.color,
        }
        return kwargs

    def fill_between(self):
        kwargs = {
            "alpha": self.alpha * 0.3,
            "color": self.color,
            "edgecolor": "none",
        }
        return kwargs


gran_styler = Styler(alpha=0.8, color="xkcd:azure", marker="s")
final_styler = Styler(alpha=0.6, color="k", marker="o")


def scatter_both(ax, x, y, G, styler, label_used="Used", label_unused="Unused"):
    ax.scatter(x[~G], y[~G], **styler.scatter(False), label=label_unused)
    ax.scatter(x[G], y[G], **styler.scatter(), label=label_used)
    return ax


def emf(tt, ax=None):
    """EMF changes as acid is added throughout a titration."""
    ttt = tt.titration
    if ax is None:
        fig, ax = plt.subplots(dpi=default.dpi)
    ax.scatter(ttt.titrant_mass * 1e3, ttt.emf, **final_styler.scatter())
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
            **gran_styler.scatter(False),
        )
        ax.scatter(
            ttt.titrant_mass[ttt.G_gran] * 1e3,
            ttt.pH_gran[ttt.G_gran],
            label="Gran — used",
            **gran_styler.scatter(),
        )
    ax.scatter(
        ttt.titrant_mass[~ttt.G_final] * 1e3,
        ttt.pH[~ttt.G_final],
        label="Final — unused",
        **final_styler.scatter(False),
    )
    ax.scatter(
        ttt.titrant_mass[ttt.G_final] * 1e3,
        ttt.pH[ttt.G_final],
        label="Final — used",
        **final_styler.scatter(),
    )
    pH_range = (
        tt.solve_kwargs["pH_range"]
        if "pH_range" in tt.solve_kwargs
        else default.pH_range
    )
    xlim = ax.get_xlim()
    ax.fill_between(
        xlim,
        *pH_range,
        **final_styler.fill_between(),
        zorder=-1,
        label="Target pH range",
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
        **gran_styler.scatter(False),
    )
    ax.scatter(
        ttt.titrant_mass[ttt.G_gran] * 1e3,
        ttt.emf0_gran[ttt.G_gran],
        label="Used",
        **gran_styler.scatter(),
    )
    ax.set_xlabel("Titrant mass / g")
    ax.set_ylabel("EMF° estimate / mV")
    ax.axhline(tt.emf0_gran, **gran_styler.plot(), label="Gran EMF°")
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
        label="Unused",
        **gran_styler.scatter(False),
    )
    ax.scatter(
        ttt.titrant_mass[ttt.G_gran] * 1e3,
        ttt.gran_estimates[ttt.G_gran],
        label="Used",
        **gran_styler.scatter(),
    )
    ax.axhline(0, c="k", lw=0.8, zorder=-1)
    x_line = np.array([-tt.gran_intercept / tt.gran_slope, ttt.titrant_mass.max()])
    y_line = x_line * tt.gran_slope + tt.gran_intercept
    ax.plot(x_line * 1e3, y_line, **gran_styler.plot(), label="Fit")
    ax.axvline(
        x_line[0] * 1e3, **gran_styler.plot(), lw=1, dashes=[3, 3], label="Intercept"
    )
    ax.set_xlabel("Titrant mass / g")
    ax.set_ylabel("Gran function")
    ax.set_title("Gran alkalinity = {:.1f} μmol/kg-sw".format(tt.alkalinity_gran))
    ax.legend()
    misc.add_credit(ax)
    return ax


def alkalinity(tt, ax=None):
    """Plot final estimates of alkalinity throughout the titration."""
    ttt = tt.titration
    if ax is None:
        fig, ax = plt.subplots(dpi=default.dpi)
    ax.scatter(
        ttt.titrant_mass[~ttt.G_final] * 1e3,
        ttt.alkalinity_estimate[~ttt.G_final] * 1e6,
        label="Unused",
        **final_styler.scatter(False),
    )
    ax.scatter(
        ttt.titrant_mass[ttt.G_final] * 1e3,
        ttt.alkalinity_estimate[ttt.G_final] * 1e6,
        label="Used",
        **final_styler.scatter(),
    )
    ax.axhline(tt.alkalinity, **final_styler.plot(), label="Final value", zorder=-1)
    ax.set_xlabel("Titrant mass / g")
    ax.set_ylabel("Alkalinity / μmol/kg-sw")
    ax.set_title(
        "Total alkalinity = ({:.1f} ± {:.1f}) μmol/kg-sw; $n$ = {}".format(
            tt.alkalinity,
            ttt.alkalinity_estimate[ttt.G_final].std() * 1e6,
            ttt.G_final.sum(),
        )
    )
    ax.legend()
    misc.add_credit(ax)
    return ax


component_stylers = {
    "alkalinity_estimate": final_styler,
    "HCO3": Styler(color="r"),
    "CO3": Styler(color="b"),
    "BOH4": Styler(color="g"),
    "PO4": Styler(color="xkcd:orange", alpha=0.8),
    "HPO4": Styler(color="xkcd:orange", alpha=0.7),
    "H3PO4": Styler(color="xkcd:orange", alpha=0.5),
    "HSO4": Styler(color="c"),
    "HF": Styler(color="y"),
    "H3SiO4": Styler(color="grey"),
    "NH3": Styler(color="g"),
    "HS": Styler(color="r"),
    "H": Styler(color="xkcd:strawberry"),
    "OH": Styler(color="xkcd:strawberry"),
    "alk_alpha": Styler(color="g"),
    "alk_beta": Styler(color="g"),
}

component_labels = {
    "alkalinity_estimate": "Total",
    "HCO3": "+[HCO$_3^-$]",
    "CO3": "+2[CO$_3^{2-}$]",
    "BOH4": "+[B(OH)$_4^-$]",
    "PO4": "+2[PO$_4^{3-}$]",
    "HPO4": "+[HPO$_4^{2-}$]",
    "H3PO4": "$-$[H$_3$PO$_4$]",
    "HSO4": "$-$[HSO$_4$]",
    "HF": "$-$[HF]",
    "H3SiO4": "+[H$_3$SiO$_4$]",
    "NH3": "+[NH$_3$]",
    "HS": "+[HS$^-$]",
    "H": "$-$[H$^+$]",
    "OH": "+[OH$^-$]",
    "alk_alpha": "$A_α$",
    "alk_beta": "$A_β$",
}


def components(tt, ax=None, log_scale=True):
    """Plot all contributors to alkalinity throughout the titration."""
    ttt = tt.titration
    if ax is None:
        fig, ax = plt.subplots(dpi=default.dpi, figsize=(6, 6))
    x_line = np.linspace(ttt.titrant_mass.min(), ttt.titrant_mass.max(), num=500)
    for var, styler in component_stylers.items():
        if var in ttt:
            if (ttt[var] > 0).all():
                if log_scale:
                    yvar = -np.log10(
                        ttt[var] * np.abs(simulate.component_multipliers[var])
                    )
                else:
                    yvar = ttt[var] * simulate.component_multipliers[var] * 1e6
                ax.plot(
                    x_line * 1e3,
                    pchip_interpolate(
                        ttt.titrant_mass.to_numpy(),
                        yvar.to_numpy(),
                        x_line,
                    ),
                    **styler.plot(),
                    label=component_labels[var] if var in component_labels else var,
                )
                scatter_both(
                    ax,
                    ttt.titrant_mass * 1e3,
                    yvar,
                    ttt.G_final,
                    styler,
                    label_used=None,
                    label_unused=None,
                )
    ax.set_xlabel("Titrant mass / g")
    if log_scale:
        ax.invert_yaxis()
        ax.set_ylabel("p(Content / mol/kg-sw)")
    else:
        ax.axhline(0, c="k", lw=0.8)
        ax.set_ylabel("Content / μmol/kg-sw")
    ax.legend(bbox_to_anchor=(1.05, 1))
    misc.add_credit(ax)
