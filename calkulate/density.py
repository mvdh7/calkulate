# Calkulate: seawater total alkalinity from titration data.
# Copyright (C) 2019-2020  Matthew Paul Humphreys  (GNU GPLv3)
"""Estimate seawater and acid densities from temperature and salinity."""
from numpy import sqrt

def sw(tempK, pSal):
    """Seawater density at 1 atm in kg/l (Millero and Poisson, 1981)."""
    # 0 < T < 40 degC & 0.5 < S < 43
    tempC = tempK - 273.15
    return (999.842594
          +   6.793952e-2 * tempC
          -   9.095290e-3 * tempC**2
          +   1.001685e-4 * tempC**3
          -   1.120083e-6 * tempC**4
          +   6.536336e-9 * tempC**5
      + (     0.824493
          -   4.0899e-3   * tempC
          +   7.6438e-5   * tempC**2
          -   8.2467e-7   * tempC**3
          +   5.3875e-9   * tempC**4 ) * pSal
      + ( -   5.72466e-3
          +   1.0227e-4   * tempC
          -   1.6546e-6   * tempC**2 ) * pSal**1.5
      +       4.8314e-4                * pSal**2   ) * 1e-3

# --- Estimate HCl titrant density --------------------------------------------
# Output in kg/l
# Uses a second-order polynomial, fit through temperature/density points:
# tempK = [288.15, 290.65, 293.15, 295.65, 298.15, 300.65,
#       303.15, 305.65, 308.15] # K
# D  = [1.025664, 1.025079, 1.024442, 1.023754, 1.023018, 1.022235,
#       1.021408, 1.020538, 1.019627] # kg/l
# These density points were calculated using E-AIM
#  [http://www.aim.env.uea.ac.uk/aim/density/density_electrolyte.php]
#  with option "rho, at the total solute mass fraction"
#  and the concentrations:
#    H+ : 0.1 mol/l
#   Na+ : 0.6 mol/l
#   Cl- : 0.7 mol/l
# This represents (approximately) a 0.1 molar HCl titrant in an NaCl solution
#  with ionic strength equal to that of seawater (cf. DAA03)
def acid(tempK):
    return (- 3.7239826839826254e-06*tempK**2
            + 1.9182242077921724e-03*tempK
            + 7.8213696227965890e-01)

# Or just at 25 degC, following Dickson et al. (2007), Chapter 5, Section 4.4:
# (note incorrect check value is returned; should be 1.02056 for (0.2, 0.5))
def acid25(mHCl, mNaCl):
    rhow25 = 0.99704 # g / cm**3
    m = mHCl + mNaCl # mol / kg-H2O
    # Eqs. (16) and (17):
    phiHCl  = 17.854 + 1.460*sqrt(m) - 0.307*m
    phiNaCl = 16.613 + 1.811*sqrt(m) + 0.094*m
    # Eqs. (14) and (15):
    mT = (36.46*mHCl + 58.44*mNaCl)/(mHCl + mNaCl)
    phimix = (mHCl*phiHCl + mNaCl*phiNaCl)/(mHCl + mNaCl)
    # Eq. (13):
    rho25 = (rhow25*(1e3 + mT*(mHCl + mNaCl)) /
        (1e3 + phimix*(mHCl + mNaCl)*rhow25))
    return rho25 # g / cm**3
