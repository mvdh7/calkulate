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
        self.get_equilibrium_constants()

    write_dat = io.write_dat

    def get_dilution_factor(self):
        """Factor for dilution of the analyte by the titrant."""
        self.mixture.dilution_factor = self.analyte.mass / (
            self.analyte.mass + self.titrant.mass
        )

    def dilute(self):
        """Calculate and apply dilution factor to total salts in the mixture.
        Assumes that salinity is not diluted, i.e. titrant and analyte have same salinity.
        """
        self.get_dilution_factor()
        # For convenience:
        mix = self.mixture
        df = mix.dilution_factor
        # User provides or assumed zero:
        mix.total_ammonia = self.analyte.total_ammonia * df
        mix.total_phosphate = self.analyte.total_phosphate * df
        mix.total_silicate = self.analyte.total_silicate * df
        mix.total_sulfide = self.analyte.total_sulfide * df
        # User provides or estimated later from salinity:
        if self.analyte.total_borate is not None:
            mix.total_borate = self.analyte.total_borate * df
        else:
            mix.total_borate = None
        if self.analyte.total_fluoride is not None:
            mix.total_fluoride = self.analyte.total_fluoride * df
        else:
            mix.total_fluoride = None
        if self.analyte.total_sulfate is not None:
            mix.total_sulfate = self.analyte.total_sulfate * df
        else:
            mix.total_sulfate = None

    def get_total_salts(self):
        """Get dict of total salt concentrations from and for PyCO2SYS."""
        conditioned, npts = pyco2.engine.condition(
            {
                "SAL": self.analyte.salinity,
                "NH3": self.mixture.total_ammonia,
                "PO4": self.mixture.total_phosphate,
                "SI": self.mixture.total_silicate,
                "H2S": self.mixture.total_sulfide,
                "K1K2CONSTANTS": self.settings.carbonic_constants,
                "BORON": self.settings.borate_ratio,
            },
        )
        # Include any internal overrides that have been provided
        totals = {}
        if self.analyte.total_borate is not None:
            totals.update({"TB": self.mixture.total_borate * 1e-6})
        if self.analyte.total_fluoride is not None:
            totals.update({"TF": self.mixture.total_fluoride * 1e-6})
        if self.analyte.total_sulfate is not None:
            totals.update({"TSO4": self.mixture.total_sulfate * 1e-6})
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

    def get_equilibrium_constants(self):
        """Get dict of equilibrium constants from and for PyCO2SYS."""
        assert hasattr(
            self.mixture, "total_salts"
        ), "You must run get_total_salts() before get_equilibrium_constants()."
        conditioned = pyco2.engine.condition(
            {
                "TEMPIN": self.mixture.temperature,
                "PRESIN": 0,
                "pHSCALEIN": 3,
                "K1K2CONSTANTS": self.settings.carbonic_constants,
                "KSO4CONSTANT": self.settings.bisulfate_constant,
                "KFCONSTANT": self.settings.fluoride_constant,
            }
        )[0]
        self.mixture.equilibrium_constants = pyco2.equilibria.assemble(
            conditioned["TEMPIN"],
            conditioned["PRESIN"],
            self.mixture.total_salts,
            conditioned["pHSCALEIN"],
            conditioned["K1K2CONSTANTS"],
            conditioned["KSO4CONSTANT"],
            conditioned["KFCONSTANT"],
        )
