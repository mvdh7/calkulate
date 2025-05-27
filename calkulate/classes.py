import pandas as pd


class Dataset(pd.DataFrame):
    """
    `calk.Dataset`
    ==============
    A `calk.Dataset` is pandas `DataFrame` with several `calk.dataset`
    functions available as methods.

    Alkalinity solving methods
    --------------------------
    calibrate
        Find the best-fit `titrant_molinity` for each sample that has a value
        for `alkalinity_certified`.
    solve
        Solve every sample for `alkalinity` when `titrant_molinity` is known.
    calkulate
        Run the `calibrate` and `solve` steps sequentially.

    Data visualisation methods
    --------------------------
    plot_titrant_molinity
        Plot the `titrant_molinity` values determined with `calibrate` through time.
    plot_alkalinity_offset
        Plot the difference between solved `alkalinity` and `alkalinity_certified`
        through time.

    Conversion methods
    ------------------
    to_Titration
        Return a `calk.Titration` object for one row in the `Dataset`.
    to_pandas
        Return a copy of the `Dataset` as a standard pandas `DataFrame`.
    """

    from .dataset import calibrate, calkulate, solve

    # get_batches = get_batches
    # get_total_salts = get_total_salts
    # to_Titration = to_Titration

    # from .plot import (
    #     alkalinity_offset as plot_alkalinity_offset,
    #     titrant_molinity as plot_titrant_molinity,
    # )

    def to_pandas(self):
        """Return a copy of the `Dataset` as a standard pandas `DataFrame`."""
        return pd.DataFrame(self)
