# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019--2021  Matthew P. Humphreys  (GNU GPLv3)
"""Calculate the density of various solutions."""

import numpy as np
from . import default


def seawater_1atm_MP81(temperature=25, salinity=35):
    """Seawater density at 1 atm in kg/l following MP81.

    Validity:
      *  0 < temperature < 40 °C
      *  0.5 < salinity < 43
    """
    if np.any(temperature < 0) or np.any(temperature > 40):
        print("Warning: some temperature values fall outside the valid range")
        print("of the MP81 density equation (0–40 °C).")
    if np.any(salinity < 0.5) or np.any(salinity > 43):
        print("Warning: some salinity values fall outside the valid range")
        print("of the MP81 density equation (0.5–43).")
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


def HCl_NaCl_25C_DSC07(
    molinity_HCl=default.molinity_HCl, molinity_NaCl=default.molinity_NaCl
):
    """Density of a mixture of HCl and NaCl at 25 °C and 1 atm following DSC07."""
    rhow25 = 0.99704  # g / cm**3
    # For convenience:
    mHCl = molinity_HCl
    mNaCl = molinity_NaCl
    molinity_total = mHCl + mNaCl  # mol / kg-H2O
    # DSC07 eqs. 16 and 17:
    phiHCl = 17.854 + 1.460 * np.sqrt(molinity_total) - 0.307 * molinity_total
    phiNaCl = 16.613 + 1.811 * np.sqrt(molinity_total) + 0.094 * molinity_total
    # DSC07 eqs. 14 and 15:
    mT = (36.46 * mHCl + 58.44 * mNaCl) / (mHCl + mNaCl)
    phi_mix = (mHCl * phiHCl + mNaCl * phiNaCl) / (mHCl + mNaCl)
    # DSC07 eq. 13:
    rho25 = (
        rhow25 * (1e3 + mT * (mHCl + mNaCl)) / (1e3 + phi_mix * (mHCl + mNaCl) * rhow25)
    )
    return rho25


def H2SO4_25C_EAIM(titrant_molinity):
    """Density of H2SO4 at 25 °C and 1 atm based on a quadratic parameterisation of data
    calculated with E-AIM: http://www.aim.env.uea.ac.uk/aim/model1/model1c.php
    """
    if (titrant_molinity < 0.05) or (titrant_molinity > 3):
        print("Warning: titrant_molinity out of parameterisation range!")
    return (
        0.99750018 + 0.06181819 * titrant_molinity - 0.00285491 * titrant_molinity ** 2
    )
