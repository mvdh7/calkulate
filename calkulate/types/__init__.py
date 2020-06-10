# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew P. Humphreys  (GNU GPLv3)
"""Classes for different types of titration."""

import numpy as np
import PyCO2SYS as pyco2
from . import components, io


class Potentiometric:
    """A potentiometrically monitored titration."""

    def __init__(self, titration_table_row, **read_dat_kwargs):
        ttr = titration_table_row  # for convenience
        self.analyte = components.Analyte()
        self.titrant = components.Titrant()
        self.mixture = components.Mixture()
        self.read_dat(ttr.fname, **read_dat_kwargs)

    # Import and export functions
    read_dat = io.read_dat
    write_dat = io.write_dat

    # def get_salts(
    #     self,
    #     salinity,
    #     ammonia=0,
    #     borate_ratio_opt=2,
    #     carbonic_constants_opt=10,
    #     phosphate=0,
    #     silicate=0,
    #     sulfide=0,
    #     totals=None,
    # ):
    #     """Get dict of total salt concentrations from and for PyCO2SYS."""
    #     conditioned = pyco2.engine.condition(
    #         {
    #             "SAL": salinity,
    #             "NH3": ammonia,
    #             "PO4": phosphate,
    #             "SI": silicate,
    #             "H2S": sulfide,
    #             "K1K2CONSTANTS": carbonic_constants_opt,
    #             "BORON": borate_ratio_opt,
    #         },
    #         npts=1,
    #     )[0]
    #     if totals is not None:
    #         totals = {k: v * 1e-6 for k, v in totals.items()}
    #         totals = pyco2.engine.condition(totals, npts=1)[0]
    #     self.analyte.total_salts = pyco2.salts.assemble(
    #         conditioned["SAL"],
    #         conditioned["SI"],
    #         conditioned["PO4"],
    #         conditioned["NH3"],
    #         conditioned["H2S"],
    #         conditioned["K1K2CONSTANTS"],
    #         conditioned["BORON"],
    #         totals=totals,
    #     )
