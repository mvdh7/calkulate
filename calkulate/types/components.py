# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew P. Humphreys  (GNU GPLv3)
"""Classes for the different components of a titration."""

import copy
import textwrap
import numpy as np
from .. import density, io


class Analyte:
    """Properties of the analyte being titrated."""

    def __init__(self, ttr, fdata):
        self.salinity = float(ttr.salinity)
        self.temperature = float(fdata["mixture_temperature"][0])
        self.mass = io.check_set(ttr, "analyte_mass", None)
        if self.mass is None:
            self.volume = io.check_set(ttr, "analyte_volume", None)
            self.density = density.seawater_atm_MP81(self.temperature, self.salinity)
            self.mass = self.volume * self.density
        self.mass *= 1e-3  # convert g to kg
        self.alkalinity_certified = io.check_set(ttr, "alkalinity_certified", None)
        if self.alkalinity_certified is not None:
            self.alkalinity_certified = float(self.alkalinity_certified)
        # User provides or assumed zero:
        self.total_ammonia = float(io.check_set(ttr, "total_ammonia", 0))
        self.total_carbonate = float(io.check_set(ttr, "total_carbonate", 0))
        self.total_phosphate = float(io.check_set(ttr, "total_phosphate", 0))
        self.total_silicate = float(io.check_set(ttr, "total_silicate", 0))
        self.total_sulfide = float(io.check_set(ttr, "total_sulfide", 0))
        # User provides or estimated later from salinity:
        self.total_borate = io.check_set(ttr, "total_borate", None)
        self.total_fluoride = io.check_set(ttr, "total_fluoride", None)
        self.total_sulfate = io.check_set(ttr, "total_sulfate", None)

    def __repr__(self):
        return textwrap.dedent(
            """\
            calkulate.types.components.Analyte()
                         mass = {:>9.3f} g
                     salinity = {:>9.3f}
                  temperature = {:>9.3f} °C
                total_ammonia = {:>9.3f} μmol/kg-sw
                 total_borate = {:>9.3f} μmol/kg-sw
              total_carbonate = {:>9.3f} μmol/kg-sw
               total_fluoride = {:>9.3f} μmol/kg-sw
              total_phosphate = {:>9.3f} μmol/kg-sw
               total_silicate = {:>9.3f} μmol/kg-sw
                total_sulfate = {:>9.3f} μmol/kg-sw
                total_sulfide = {:>9.3f} μmol/kg-sw\
            """.format(
                self.mass * 1e3,
                self.salinity,
                self.temperature,
                self.total_ammonia,
                self.total_borate,
                self.total_carbonate,
                self.total_fluoride,
                self.total_phosphate,
                self.total_silicate,
                self.total_sulfate,
                self.total_sulfide,
            )
        )


class Titrant:
    """Properties of the titrant added to the analyte."""

    def __init__(self, ttr, fdata):
        amount_unit = io.check_set(ttr, "titrant_amount_unit", "ml")
        assert amount_unit in ["ml", "g"]
        if amount_unit == "ml":
            self.volume = fdata["titrant_amount"]
            self.density = density.HCl_NaCl_25C_DSC07()
            self.mass = self.volume * self.density * 1e-3
        elif amount_unit == "g":
            self.mass = fdata["titrant_amount"] * 1e-3
        self.increments = np.size(self.mass)
        self.concentration = io.check_set(ttr, "titrant_concentration", None)
        self.molinity = io.check_set(ttr, "titrant_molinity", None)
        if self.concentration is not None and self.molinity is None:
            self.set_molinity()

    def set_molinity(self):
        self.molinity = self.concentration / self.density

    def subset(self, use_points=None):
        """Return a subset of the Titrant."""
        subtitrant = copy.deepcopy(self)
        subtitrant.mass = self.mass[use_points].ravel()
        if hasattr(self, "volume"):
            subtitrant.volume = self.volume[use_points].ravel()
        return subtitrant

    # def __repr__(self):
    #     return textwrap.dedent(
    #         """\
    #         calkulate.types.components.Titrant
    #                molinity = {:>5.3f} mol/kg
    #                    mass = from {:>5.3f} to {:>5.3f} g
    #              increments = {}\
    #         """.format(
    #             self.molinity,
    #             np.min(self.mass) * 1e3,
    #             np.max(self.mass) * 1e3,
    #             self.increments,
    #         )
    #     )


class Mixture:
    """Properties of the titrant-analyte mixture."""

    def __init__(self, ttr, fdata, measurement_type="EMF"):
        if measurement_type == "EMF":
            self.emf = fdata["mixture_measurement"]
        elif measurement_type == "pH":
            self.pH = fdata["mixture_measurement"]
        self.temperature = fdata["mixture_temperature"]
        if "temperature_override" in ttr:
            if ~np.isnan(ttr.temperature_override):
                self.temperature[:] = float(ttr.temperature_override)
        # Equilibrium constant overrides
        self.k_ammonia = io.check_set(ttr, "k_ammonia", None)
        self.k_boric = io.check_set(ttr, "k_boric", None)
        self.k_carbonic_1 = io.check_set(ttr, "k_carbonic_1", None)
        self.k_carbonic_2 = io.check_set(ttr, "k_carbonic_2", None)
        self.k_hydrofluoric = io.check_set(ttr, "k_hydrofluoric", None)
        self.k_phosphoric_1 = io.check_set(ttr, "k_phosphoric_1", None)
        self.k_phosphoric_2 = io.check_set(ttr, "k_phosphoric_2", None)
        self.k_phosphoric_3 = io.check_set(ttr, "k_phosphoric_3", None)
        self.k_orthosilicic = io.check_set(ttr, "k_orthosilicic", None)
        self.k_hydrosulfuric_1 = io.check_set(ttr, "k_hydrosulfuric_1", None)
        self.k_sulfuric_2 = io.check_set(ttr, "k_sulfuric_2", None)
        self.k_water = io.check_set(ttr, "k_water", None)
        self.equilibrium_constants = None

    def subset(self, use_points=None):
        """Return a subset of the Mixture."""
        submixture = copy.deepcopy(self)
        for k, v in submixture.__dict__.items():
            if k in ["total_salts", "equilibrium_constants"]:
                submixture.__dict__[k] = {
                    m: w[use_points].ravel() for m, w in submixture.__dict__[k].items()
                }
            elif k.startswith("k_"):
                submixture.__dict__[k] = v
            else:
                submixture.__dict__[k] = v[use_points].ravel()
        return submixture

    def _overwrite_equilibrium_constant(self, equilibrium_constant, pyco2_constant):
        overwrite_value = self.__dict__[equilibrium_constant]
        if overwrite_value is not None:
            self.equilibrium_constants[pyco2_constant][:] = overwrite_value

    def __repr__(self):
        return textwrap.dedent(
            """\
            calkulate.types.components.Mixture
              Attributes: {}\
            """.format(
                list(self.__dict__.keys())
            )
        )


class Settings:
    """Settings for solving the titration."""

    def __init__(self, ttr, fdata):
        self.bisulfate_constant = int(io.check_set(ttr, "bisulfate_constant", 1))
        self.borate_ratio = int(io.check_set(ttr, "borate_ratio", 2))
        self.carbonic_constants = int(io.check_set(ttr, "carbonic_constants", 10))
        self.fluoride_constant = int(io.check_set(ttr, "fluoride_constant", 1))

    def __repr__(self):
        return textwrap.dedent(
            """\
            calkulate.types.components.Settings
              bisulfate_constant = {:>2}
                    borate_ratio = {:>2}
              carbonic_constants = {:>2}
               fluoride_constant = {:>2}\
            """.format(
                self.bisulfate_constant,
                self.borate_ratio,
                self.carbonic_constants,
                self.fluoride_constant,
            )
        )
