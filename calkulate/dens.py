from .const import Tzero

# Seawater density in kg/l   
# Inputs: T (temperature / deg C), S (practical salinity)
# Valid ranges: 0 < T < 40 degC & 0.5 < S < 43
# Source: Millero & Poisson (1981) DSR-A 28(6).
#          doi:10.1016/0198-0149(81)90122-9

def sw(Tk,S):
    
    Tc = Tk - Tzero
    
    return (999.842594                       \
          +   6.793952e-2 * Tc               \
          -   9.095290e-3 * Tc**2            \
          +   1.001685e-4 * Tc**3            \
          -   1.120083e-6 * Tc**4            \
          +   6.536336e-9 * Tc**5            \
      + (     0.824493                       \
          -   4.0899e-3   * Tc               \
          +   7.6438e-5   * Tc**2            \
          -   8.2467e-7   * Tc**3            \
          +   5.3875e-9   * Tc**4 ) * S      \
      + ( -   5.72466e-3                     \
          +   1.0227e-4   * Tc               \
          -   1.6546e-6   * Tc**2 ) * S**1.5 \
      +       4.8314e-4             * S**2   ) * 1e-3

# =============================================================================
#
# --- Estimate HCl titrant density --------------------------------------------
#
# Output in kg/l
#      
# Uses a second-order polynomial, fit through temperature/density points:
#    
# Tk = [288.15, 290.65, 293.15, 295.65, 298.15, 300.65, 
#       303.15, 305.65, 308.15] # K
# D  = [1.025664, 1.025079, 1.024442, 1.023754, 1.023018, 1.022235,
#       1.021408, 1.020538, 1.019627] # kg/l
#
# These density points were calculated using E-AIM
#  [http://www.aim.env.uea.ac.uk/aim/density/density_electrolyte.php]
#  with option "rho, at the total solute mass fraction"
#  and the concentrations:
#    H+ : 0.1 mol/l
#   Na+ : 0.6 mol/l
#   Cl- : 0.7 mol/l
# This represents (approximately) a 0.1 molar HCl titrant in an NaCl solution
#  with ionic strength equal to that of seawater (cf. DAA03)
    
def acid(Tk):

    return - 3.7239826839826254e-06 * Tk**2 \
           + 1.9182242077921724e-03 * Tk    \
           + 7.8213696227965890e-01
      