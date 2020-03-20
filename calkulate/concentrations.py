# Calkulate: seawater total alkalinity from titration data.
# Copyright (C) 2019-2020  Matthew Paul Humphreys  (GNU GPLv3)
"""Estimate total concentrations from practical salinity."""
from PyCO2SYS import assemble

def concTotals(pSal, totalCarbonate=0, totalPhosphate=0, totalSilicate=0,
        totalAmmonia=0, totalH2Sulfide=0, WhichKs=10, WhoseTB=2):
    """Assemble a dict of sample concentrations in mol/kg-sw using PyCO2SYS."""
    aconcs = assemble.inputs({
        'SAL': pSal, 'WhichKs': WhichKs, 'WhoseTB': WhoseTB,})[0].values()
    totalBorate, totalFluoride, totalSulfate = \
        assemble.concentrations(*aconcs)
    return {
        'C': totalCarbonate,
        'B': totalBorate,
        'S': totalSulfate,
        'F': totalFluoride,
        'P': totalPhosphate,
        'Si': totalSilicate,
        'NH3': totalAmmonia,
        'H2S': totalH2Sulfide,
        }
