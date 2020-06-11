# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew P. Humphreys  (GNU GPLv3)
"""Import and export titration data."""

import numpy as np


def read_dat(
    fname,
    acid_volume_col=0,
    emf_col=1,
    temperature_col=2,
    delimiter="\t",
    skip_header=2,
    **kwargs
):
    """Import potentiometric titration data from a text file."""
    data = np.genfromtxt(fname, delimiter=delimiter, skip_header=skip_header, **kwargs)
    titrant_volume = data[:, acid_volume_col]
    mixture_emf = data[:, emf_col]
    mixture_temperature = data[:, temperature_col]
    return {
        "titrant_volume": titrant_volume,
        "mixture_emf": mixture_emf,
        "mixture_temperature": mixture_temperature,
    }


def write_dat(
    titration,
    fname,
    line0="",
    line1="titrant_volume\temf\ttemperature",
    mode="x",
    volume_fmt=".5f",
    emf_fmt=".5f",
    temperature_fmt=".3f",
    **kwargs
):
    """Write potentiometric titration data to a text file."""
    with open(fname, mode=mode, **kwargs) as f:
        f.write("{}\n{}\n".format(line0, line1))
        for i in range(len(titration.titrant.volume)):
            f.write(
                (
                    "{:"
                    + volume_fmt
                    + "}\t{:"
                    + emf_fmt
                    + "}\t{:"
                    + temperature_fmt
                    + "}\n"
                ).format(
                    titration.titrant.volume[i],
                    titration.mixture.emf[i],
                    titration.mixture.temperature[i],
                )
            )


def check_set(collection, field, default):
    """Check if a field is in a collection and that it is not a NaN.
    Return the field if so, or the default if not.
    """
    if field in collection:
        if ~np.isnan(collection[field]):
            return collection[field]
        else:
            return default
    else:
        return default
