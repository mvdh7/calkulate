# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew P. Humphreys  (GNU GPLv3)
"""Classes for different types of titration."""

import numpy as np
from . import components, io


class Potentiometric:
    """A potentiometrically monitored titration."""

    def __init__(self, fname=None, **read_dat_kwargs):
        self.analyte = components.Analyte()
        self.titrant = components.Titrant()
        self.mixture = components.Mixture()
        if fname is not None:
            self.read_dat(fname, **read_dat_kwargs)

    # Import and export functions
    read_dat = io.read_dat
    write_dat = io.write_dat
