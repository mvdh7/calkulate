# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew P. Humphreys  (GNU GPLv3)
import numpy as np
from . import constants


def gran_estimator(
    mixture_mass, emf, temperature, use_points=None,
):
    """Simple Gran-plot estimator following DAA03 eq. 10."""
    if use_points is None:
        G = np.full(np.size(emf), True)
    else:
        G = use_points
    temperature_K = temperature + constants.absolute_zero
    return mixture_mass[G] * np.exp(
        emf[G] * constants.faraday / (constants.ideal_gas * temperature_K[G])
    )
