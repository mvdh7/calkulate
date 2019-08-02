# Calkulate: seawater total alkalinity from titration data.
# Copyright (C) 2019  Matthew Paul Humphreys  (GNU GPLv3)
"""Estimate total concentrations from practical salinity."""
from .constants import pSal2cl, ramB, ramF

def XT(pSal, CT=0, PT=0, SiT=0):
    """Assemble a dict of sample concentrations in mol/kg-sw."""
    return {
        'C': CT,
        'B': BT(pSal),
        'S': ST(pSal),
        'F': FT(pSal),
        'P': PT,
        'Si': SiT,
    } # all in mol/kg-sw

def BT(pSal):
    """Estimate total borate from practical salinity in mol/kg-sw [LKB10]."""
    return pSal*0.1336e-3/ramB

def FT(pSal):
    """Estimate total fluoride from practical salinity in mol/kg-sw [W71]."""
    return pSal*6.75e-5/(ramF*pSal2cl)

def ST(pSal):
    """Estimate total sulfate from practical salinity in mol/kg-sw [???]."""
    return (0.14/96.061)*(pSal/pSal2cl)
