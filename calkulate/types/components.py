# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew P. Humphreys  (GNU GPLv3)
"""Classes for the different components of a titration."""

import textwrap
import numpy as np
from .. import density, io


class Analyte:
    """Properties of the analyte being titrated."""

    def __init__(self, ttr, fdata):
        self.salinity = float(ttr.salinity)
        self.temperature = float(fdata["mixture_temperature"][0])
        self.density = density.seawater_atm_MP81(self.temperature, self.salinity)
        self.volume = float(ttr.analyte_volume)
        self.mass = self.volume * self.density
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
                      density = {:>9.3f} kg/l
                         mass = {:>9.3f} g
                     salinity = {:>9.3f}
                  temperature = {:>9.3f} °C
                       volume = {:>9.3f} ml
                total_ammonia = {:>9.3f} μmol/kg-sw
                 total_borate = {:>9.3f} μmol/kg-sw
              total_carbonate = {:>9.3f} μmol/kg-sw
               total_fluoride = {:>9.3f} μmol/kg-sw
              total_phosphate = {:>9.3f} μmol/kg-sw
               total_silicate = {:>9.3f} μmol/kg-sw
                total_sulfate = {:>9.3f} μmol/kg-sw
                total_sulfide = {:>9.3f} μmol/kg-sw\
            """.format(
                self.density,
                self.mass,
                self.salinity,
                self.temperature,
                self.volume,
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
        self.volume = fdata["titrant_volume"]
        self.density = density.HCl_NaCl_25C_DSC07()
        self.mass = self.volume * self.density
        self.concentration = io.check_set(ttr, "titrant_concentration", None)
        if self.concentration is not None:
            self.set_molality()
        self.increments = np.size(self.mass)

    def set_molality(self):
        self.molality = self.concentration / self.density

    def __repr__(self):
        return textwrap.dedent(
            """\
            calkulate.types.components.Titrant
              concentration = {:>5.3f} mol/l
                    density = {:>5.3f} kg/l
                 increments = {}
                       mass = from {:>5.3f} to {:>5.3f} g
                   molality = {:>5.3f} mol/kg
                     volume = from {:>5.3f} to {:>5.3f} ml\
            """.format(
                self.concentration,
                self.density,
                self.increments,
                np.min(self.mass),
                np.max(self.mass),
                self.molality,
                np.min(self.volume),
                np.max(self.volume),
            )
        )


class Mixture:
    """Properties of the titrant-analyte mixture."""

    def __init__(self, ttr, fdata):
        self.emf = fdata["mixture_emf"]
        self.temperature = fdata["mixture_temperature"]

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
