# Calkulate: seawater total alkalinity from titration data.
# Copyright (C) 2019  Matthew Paul Humphreys  (GNU GPLv3)
"""Estimate total concentrations from practical salinity."""
from .constants import pSal2cl, ramB, ramF

def concTotals(pSal, totalCarbonate=0, totalPhosphate=0, totalSilicate=0):
    """Assemble a dict of sample concentrations in mol/kg-sw."""
    return {
        'C': totalCarbonate,
        'B': totalBorate_LKB10(pSal),
        'S': totalSulfate_MR66(pSal),
        'F': totalFluoride_W71(pSal),
        'P': totalPhosphate,
        'Si': totalSilicate,
    }

def totalBorate_LKB10(pSal):
    """Estimate total borate from practical salinity in mol/kg-sw [LKB10]."""
    return pSal*0.1336e-3/ramB

def totalBorate_U74(pSal):
    """Estimate total borate from practical salinity in mol/kg-sw [U74]."""
    return pSal*0.232e-3/(ramB*pSal2cl)

def totalFluoride_R65(pSal):
    """Estimate total fluoride from practical salinity in mol/kg-sw [R65]."""
    return pSal*6.7e-5/(ramF*pSal2cl)
    
def totalFluoride_W71(pSal):
    """Estimate total fluoride from practical salinity in mol/kg-sw [W71]."""
    return pSal*6.75e-5/(ramF*pSal2cl)

def totalSulfate_MR66(pSal):
    """Estimate total sulfate from practical salinity in mol/kg-sw [MR66]."""
    return (0.14/96.061)*(pSal/pSal2cl)
