# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019  Matthew Paul Humphreys  (GNU GPLv3)

from .constants import psal2Cl, RMM_B, RMM_F

def XT(psal, CT=0, PT=0, SiT=0):
    """Assemble a dict of concentrations in mol/kg-sw."""
    return {
        'C': CT,
        'B': BT(psal),
        'S': ST(psal),
        'F': FT(psal),
        'P': PT,
        'Si': SiT,
    } # all in mol/kg-sw

def BT(psal):
    """Estimate total borate from practical salinity in mol/kg-sw [LKB10]."""
    return psal * 0.1336e-3 / RMM_B

def FT(psal):
    """Estimate total fluoride from practical salinity in mol/kg-sw [W71]."""
    return psal * 6.75e-5 / (RMM_F * psal2Cl)

def ST(psal):
    """Estimate total sulfate from practical salinity in mol/kg-sw [???]."""
    return (0.14 / 96.061) * (psal / psal2Cl)
