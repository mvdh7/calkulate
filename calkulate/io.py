# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew P. Humphreys  (GNU GPLv3)
"""Import and export data."""
import numpy as np


def read_dat(
    fname,
    acid_volume_col=0,
    delimiter="\t",
    emf_col=1,
    skip_header=2,
    temperature_col=2,
    **kwargs
):
    """Import potentiometric titration data from a text file."""
    data = np.genfromtxt(fname, delimiter=delimiter, skip_header=skip_header, **kwargs)
    acid_volume = data[:, acid_volume_col]
    emf = data[:, emf_col]
    temperature = data[:, temperature_col]
    return acid_volume, emf, temperature


def write_dat(
    fname, acid_volume, emf, temperature, line0="", line1="", mode="x", **kwargs
):
    """Write potentiometric titration data to a text file."""
    with open(fname, mode=mode, **kwargs) as f:
        f.write("{}\n{}\n".format(line0, line1))
        for i in range(len(acid_volume)):
            f.write(
                "{:.5f}\t{:.5f}\t{:.3f}\n".format(
                    acid_volume[i], emf[i], temperature[i]
                )
            )
