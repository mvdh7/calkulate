# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew P. Humphreys  (GNU GPLv3)
"""Classes for different types of titration."""

import numpy as np
import PyCO2SYS as pyco2
from .. import io
from . import components


class Potentiometric:
    """A potentiometrically monitored titration."""

    def __init__(self, titration_table_row, **read_dat_kwargs):
        ttr = titration_table_row  # for convenience
        fdata = io.read_dat(ttr.fname, **read_dat_kwargs)
        self.fname = ttr.fname
        self.settings = components.Settings(ttr, fdata)
        self.analyte = components.Analyte(ttr, fdata)
        self.titrant = components.Titrant(ttr, fdata)
        self.mixture = components.Mixture(ttr, fdata)
        self.dilute()
        self.get_total_salts()

    write_dat = io.write_dat

    def get_dilution_factor(self):
        """Factor for dilution of the analyte by the titrant."""
        self.mixture.dilution_factor = self.analyte.mass / (
            self.analyte.mass + self.titrant.mass
        )

    def dilute(self):
        """Calculate and apply dilution factor to total salts in the mixture."""
        self.get_dilution_factor()
        # For convenience:
        mix = self.mixture
        df = mix.dilution_factor
        # User provides or assumed zero:
        mix.salinity = self.analyte.salinity * df
        mix.ammonia = self.analyte.ammonia * df
        mix.phosphate = self.analyte.phosphate * df
        mix.silicate = self.analyte.silicate * df
        mix.sulfide = self.analyte.sulfide * df
        # User provides or estimated later from salinity:
        if self.analyte.borate is not None:
            mix.borate = self.analyte.borate * df
        else:
            mix.borate = None
        if self.analyte.fluoride is not None:
            mix.fluoride = self.analyte.fluoride * df
        else:
            mix.fluoride = None
        if self.analyte.sulfate is not None:
            mix.sulfate = self.analyte.sulfate * df
        else:
            mix.sulfate = None

    def get_total_salts(self):
        """Get dict of total salt concentrations from and for PyCO2SYS."""
        conditioned, npts = pyco2.engine.condition(
            {
                "SAL": self.mixture.salinity,
                "NH3": self.mixture.ammonia,
                "PO4": self.mixture.phosphate,
                "SI": self.mixture.silicate,
                "H2S": self.mixture.sulfide,
                "K1K2CONSTANTS": self.settings.carbonic_constants,
                "BORON": self.settings.borate_ratio,
            },
        )
        # Include any internal overrides that have been provided
        totals = {}
        if self.analyte.borate is not None:
            totals.update({"TB": self.mixture.borate * 1e-6})
        if self.analyte.fluoride is not None:
            totals.update({"TF": self.mixture.fluoride * 1e-6})
        if self.analyte.sulfate is not None:
            totals.update({"TSO4": self.mixture.sulfate * 1e-6})
        if len(totals) > 0:
            totals = pyco2.engine.condition(totals, npts=npts)[0]
        else:
            totals = None
        self.mixture.total_salts = pyco2.salts.assemble(
            conditioned["SAL"],
            conditioned["SI"],
            conditioned["PO4"],
            conditioned["NH3"],
            conditioned["H2S"],
            conditioned["K1K2CONSTANTS"],
            conditioned["BORON"],
            totals=totals,
        )
        