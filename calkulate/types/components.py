# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew P. Humphreys  (GNU GPLv3)
"""Classes for the different components of a titration."""

import numpy as np
from .. import density


class Analyte:
    """Properties of the analyte being titrated."""

    pass

    # def __init__(self, salinity=None, temperature=None, volume=None):
    #     if salinity is not None:
    #         assert np.size(salinity) == 1, "Analyte salinity must be scalar."
    #         assert salinity >= 0, "Analyte salinity must be positive."
    #     self.salinity = salinity
    #     if temperature is not None:
    #         assert np.size(temperature) == 1, "Analyte temperature must be scalar."
    #     self.temperature = temperature
    #     if volume is not None:
    #         assert np.size(volume) == 1, "Analyte volume must be scalar."
    #         assert volume >= 0, "Analyte volume must be positive."
    #     self.volume = volume
    #     if salinity is not None and temperature is not None:
    #         self.get_density()

    # def get_density(self):
    #     """Calculate analyte density."""
    #     self.density = density.seawater_atm_MP81(self.temperature, self.salinity)


class Titrant:
    """Properties of the titrant added to the analyte."""

    pass


class Mixture:
    """Properties of the titrant-analyte mixture."""

    pass
