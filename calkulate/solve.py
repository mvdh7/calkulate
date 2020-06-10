# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew P. Humphreys  (GNU GPLv3)
"""Solve titration data for total alkalinity."""

import numpy as np
from scipy.stats import linregress
from . import constants


def gran_estimator(titration):
    """Simple Gran-plot estimator function following DAA03 eq. 10."""
    tt = titration  # for convenience
    temperature_K = tt.mixture.temperature + constants.absolute_zero
    return tt.mixture.mass * np.exp(
        tt.mixture.emf * constants.faraday / (constants.ideal_gas * temperature_K)
    )


def gran_guess_alkalinity(titration):
    """Simple Gran-plot first-guess of alkalinity."""
    tt = titration  # for convenience
    gradient, intercept_y, _, _, _ = linregress(tt.titrant.mass, gran_estimator(tt))
    intercept_x = -intercept_y / gradient
    alkalinity_guess = intercept_x * tt.titrant.concentration / tt.analyte.mass
    return alkalinity_guess


def gran_guess_emf0(titration, HF=0, HSO4=0):
    """Simple Gran-plot first-guess of EMF0 following DAA03 eq. 11."""
    tt = titration  # for convenience
    alkalinity_guess = gran_guess_alkalinity(tt)
    return tt.mixture.emf - (constants.ideal_gas * tempK / constants.faraday) * np.log(
        (
            (tt.titrant.mass * tt.titrant.molality - tt.analyte.mass * alkalinity_guess)
            - tt.analyte.mass * (HF + HSO4)
        )
        / tt.mixture.mass
    )
