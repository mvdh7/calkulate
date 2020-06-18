# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew P. Humphreys  (GNU GPLv3)
"""Classes for different types of titration."""

import textwrap
import numpy as np
import PyCO2SYS as pyco2
from .. import convert, io, solve
from . import components


class Titration:
    """A titration dataset."""

    def __init__(self, titration_table_row, **read_dat_kwargs):
        ttr = titration_table_row  # for convenience
        self.file_path = io.check_set(ttr, "file_path", "")
        self.file_name = ttr.file_name
        fname = self.file_path + self.file_name
        fdata = io.read_dat(fname, **read_dat_kwargs)
        self.fname = fname
        self.titration_table_index = ttr.name
        self.measurement_type = io.check_set(ttr, "measurement_type", "EMF")
        self.solver = convert.measurementType_to_solver(self.measurement_type)
        self.settings = components.Settings(ttr, fdata)
        self.analyte = components.Analyte(ttr, fdata)
        self.titrant = components.Titrant(ttr, fdata)
        self.mixture = components.Mixture(
            ttr, fdata, measurement_type=self.measurement_type
        )
        self.get_mixture_mass()
        self.dilute()
        self.get_total_salts()
        self.get_equilibrium_constants()

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
        # Overwrite if user supplied values
        self.mixture._overwrite_equilibrium_constant("k_ammonia", "KNH3")
        self.mixture._overwrite_equilibrium_constant("k_boric", "KB")
        self.mixture._overwrite_equilibrium_constant("k_carbonic_1", "K1")
        self.mixture._overwrite_equilibrium_constant("k_carbonic_2", "K2")
        self.mixture._overwrite_equilibrium_constant("k_hydrofluoric", "KF")
        self.mixture._overwrite_equilibrium_constant("k_phosphoric_1", "KP1")
        self.mixture._overwrite_equilibrium_constant("k_phosphoric_2", "KP2")
        self.mixture._overwrite_equilibrium_constant("k_phosphoric_3", "KP3")
        self.mixture._overwrite_equilibrium_constant("k_orthosilicic", "KSi")
        self.mixture._overwrite_equilibrium_constant("k_hydrosulfuric_1", "KH2S")
        self.mixture._overwrite_equilibrium_constant("k_sulfuric_2", "KSO4")
        self.mixture._overwrite_equilibrium_constant("k_water", "KW")

    def get_gran_guesses(self):
        """Get initial Gran-plot guesses of alkalinity and EMF0."""
        (
            alkalinity_guess,
            self.analyte.emf0_guess,
            self.analyte.use_points_guess,
        ) = solve.gran_guesses(self)
        self.analyte.alkalinity_guess = alkalinity_guess * 1e6

    def calibrate(self, **kwargs):
        """Calibrate the titrant molinity."""
        assert (
            self.analyte.alkalinity_certified is not None
        ), "You can only calibrate if the certified alkalinity has been set."
        self.calibrated = solve.calibrate(self, solver=self.solver, **kwargs)
        self.titrant.molinity_calibrated = self.calibrated["x"][0]

    def set_own_titrant_molinity(self):
        """Set the titrant molinity to the value found by calibrating this sample."""
        self.titrant.molinity = self.calibrated["x"][0]

    def get_initial_pH(self):
        """Calculate analyte pH from EMF before any titrant addition."""
        self.analyte.pH_temperature = self.mixture.temperature[0]
        self.analyte.pH = convert.emf_to_pH(
            self.mixture.emf[0], self.analyte.emf0, self.analyte.pH_temperature
        )

    def solve(self, **kwargs):
        """Solve for total alkalinity and EMF0."""
        assert (
            self.titrant.molinity is not None
        ), "You can only solve if the titrant molinity has been set."
        self.solved = self.solver(self, **kwargs)
        self.analyte.alkalinity = self.solved["x"][0] * 1e6
        if self.measurement_type == "EMF":
            self.analyte.emf0 = self.solved["x"][1]
            self.get_initial_pH()

    def __repr__(self):
        return textwrap.dedent(
            """\
            calkulate.types.Titration
                          fname = '{}'
               Other attributes = ['analyte', 'mixture', 'settings', 'titrant']
                        Methods = write_dat()\
            """.format(
                self.fname
            )
        )
