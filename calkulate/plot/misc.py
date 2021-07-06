from .. import meta


def add_credit(ax):
    """Add Calkulate credit to figures."""
    ax.text(
        1.005,
        0,
        "Calkulate v{}".format(meta.__version__),
        alpha=0.2,
        c="k",
        ha="left",
        va="bottom",
        rotation=-90,
        transform=ax.transAxes,
    )
