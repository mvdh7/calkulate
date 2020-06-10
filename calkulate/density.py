# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew P. Humphreys  (GNU GPLv3)
"""The densities of various solutions."""


def seawater_atm_MP81(temperature, salinity):
    """Seawater density at 1 atm in kg/l following MP81."""
    # Validity: 0 < T < 40 degC & 0.5 < S < 43
    return (
        999.842594
        + 6.793952e-2 * temperature
        - 9.095290e-3 * temperature ** 2
        + 1.001685e-4 * temperature ** 3
        - 1.120083e-6 * temperature ** 4
        + 6.536336e-9 * temperature ** 5
        + (
            0.824493
            - 4.0899e-3 * temperature
            + 7.6438e-5 * temperature ** 2
            - 8.2467e-7 * temperature ** 3
            + 5.3875e-9 * temperature ** 4
        )
        * salinity
        + (-5.72466e-3 + 1.0227e-4 * temperature - 1.6546e-6 * temperature ** 2)
        * salinity ** 1.5
        + 4.8314e-4 * salinity ** 2
    ) * 1e-3
