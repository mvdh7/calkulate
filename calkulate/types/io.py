# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew P. Humphreys  (GNU GPLv3)
"""Import and export titration data."""

import numpy as np


def read_dat(
    self,
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
    self.titrant.volume = data[:, acid_volume_col]
    self.mixture.emf = data[:, emf_col]
    self.mixture.temperature = data[:, temperature_col]


def write_dat(self, fname, line0="", line1="", mode="x", **kwargs):
    """Write potentiometric titration data to a text file."""
    with open(fname, mode=mode, **kwargs) as f:
        f.write("{}\n{}\n".format(line0, line1))
        for i in range(len(self.titrant.volume)):
            f.write(
                "{:.5f}\t{:.5f}\t{:.3f}\n".format(
                    self.titrant.volume[i],
                    self.mixture.emf[i],
                    self.mixture.temperature[i],
                )
            )
