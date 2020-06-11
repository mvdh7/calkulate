# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew P. Humphreys  (GNU GPLv3)
"""Classes for different types of titration."""

import textwrap
import numpy as np
import PyCO2SYS as pyco2
from .. import convert, io, solve
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
        self.get_mixture_mass()
        self.dilute()
        self.get_total_salts()
        self.get_equilibrium_constants()
        self.get_gran_guesses()
        self.solve()

    write_dat = io.write_dat

    def get_mixture_mass(self):
        """Calculate the total mass of the titrant-analyte mixture."""
        self.mixture.mass = self.analyte.mass + self.titrant.mass

    def get_dilution_factor(self):
        """Factor for dilution of the analyte by the titrant."""
        self.mixture.dilution_factor = convert.dilution_factor(
            self.analyte.mass, self.mixture.mass
        )

    def dilute(self):
        """Calculate and apply dilution factor to total salts in the mixture.
        Assumes that salinity is not diluted, i.e. titrant and analyte have same salinity.
        """
        self.get_dilution_factor()
        df = self.mixture.dilution_factor  # for convenience
        # User provides or assumed zero:
        self.mixture.total_ammonia = self.analyte.total_ammonia * df
        self.mixture.total_carbonate = self.analyte.total_carbonate * df
        self.mixture.total_phosphate = self.analyte.total_phosphate * df
        self.mixture.total_silicate = self.analyte.total_silicate * df
        self.mixture.total_sulfide = self.analyte.total_sulfide * df
        self.mixture.salinity_analyte = self.analyte.salinity * df
        # User provides or estimated later from salinity:
        if self.analyte.total_borate is not None:
            self.mixture.total_borate = self.analyte.total_borate * df
        else:
            self.mixture.total_borate = None
        if self.analyte.total_fluoride is not None:
            self.mixture.total_fluoride = self.analyte.total_fluoride * df
        else:
            self.mixture.total_fluoride = None
        if self.analyte.total_sulfate is not None:
            self.mixture.total_sulfate = self.analyte.total_sulfate * df
        else:
            self.mixture.total_sulfate = None

    def get_total_salts(self):
        """Get dict of total salt concentrations from and for PyCO2SYS."""
        conditioned, npts = pyco2.engine.condition(
            {
                "SAL": np.insert(
                    self.mixture.salinity_analyte, 0, self.analyte.salinity
                ),
                "NH3": np.insert(
                    self.mixture.total_ammonia, 0, self.analyte.total_ammonia
                ),
                "PO4": np.insert(
                    self.mixture.total_phosphate, 0, self.analyte.total_phosphate
                ),
                "SI": np.insert(
                    self.mixture.total_silicate, 0, self.analyte.total_silicate
                ),
                "H2S": np.insert(
                    self.mixture.total_sulfide, 0, self.analyte.total_sulfide
                ),
                "K1K2CONSTANTS": self.settings.carbonic_constants,
                "BORON": self.settings.borate_ratio,
            },
        )
        # Include any internal overrides that have been provided
        totals = {}
        if self.analyte.total_borate is not None:
            totals.update(
                {
                    "TB": np.insert(
                        self.mixture.total_borate, 0, self.analyte.total_borate
                    )
                    * 1e-6
                }
            )
        if self.analyte.total_fluoride is not None:
            totals.update(
                {
                    "TF": np.insert(
                        self.mixture.total_fluoride, 0, self.analyte.total_fluoride
                    )
                    * 1e-6
                }
            )
        if self.analyte.total_sulfate is not None:
            totals.update(
                {
                    "TSO4": np.insert(
                        self.mixture.total_sulfate, 0, self.analyte.total_sulfate
                    )
                    * 1e-6
                }
            )
        if len(totals) > 0:
            totals = pyco2.engine.condition(totals, npts=npts)[0]
        else:
            totals = None
        all_total_salts = pyco2.salts.assemble(
            conditioned["SAL"],
            conditioned["SI"],
            conditioned["PO4"],
            conditioned["NH3"],
            conditioned["H2S"],
            conditioned["K1K2CONSTANTS"],
            conditioned["BORON"],
            totals=totals,
        )
        self.mixture.total_salts = {k: v[1:] for k, v in all_total_salts.items()}
        self.analyte.total_borate = all_total_salts["TB"][0] * 1e6
        self.mixture.total_borate = all_total_salts["TB"][1:] * 1e6
        self.analyte.total_fluoride = all_total_salts["TF"][0] * 1e6
        self.mixture.total_fluoride = all_total_salts["TF"][1:] * 1e6
        self.analyte.total_sulfate = all_total_salts["TSO4"][0] * 1e6
        self.mixture.total_sulfate = all_total_salts["TSO4"][1:] * 1e6

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

    def get_gran_guesses(self):
        """Get initial Gran-plot guesses of alkalinity and EMF0."""
        (
            alkalinity_guess,
            self.analyte.emf0_guess,
            self.analyte.use_points_guess,
        ) = solve.gran_guesses(self)
        self.analyte.alkalinity_guess = alkalinity_guess * 1e6

    def solve(self):
        """Solve for total alkalinity and EMF0."""
        self.solved = solve.complete(self)

    def __repr__(self):
        return textwrap.dedent(
            """\
            calkulate.types.Potentiometric
                          fname = '{}'
               Other attributes = ['analyte', 'mixture', 'settings', 'titrant']
                        Methods = write_dat()\
            """.format(
                self.fname
            )
        )
