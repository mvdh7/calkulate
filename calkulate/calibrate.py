# Calkulate: seawater total alkalinity from titration data.
# Copyright (C) 2019  Matthew Paul Humphreys  (GNU GPLv3)
"""Calibrate acid titrant concentration using certified reference material
measurements.
"""
from scipy.optimize import least_squares as olsq
from . import solve

def complete(massAcid, emf, tempK, massSample, alkCert, concTotals,
        eqConstants):
    """Calibrate acid concentration using the complete calculation method."""
    return olsq(lambda concAcid: solve.complete(massAcid, emf, tempK,
        massSample, concAcid, concTotals, eqConstants)['x'][0] - alkCert, 0.1,
        method='lm')

def DAA03(massAcid, emf, tempK, massSample, alkCert, concTotals, eqConstants):
    """Calibrate acid concentration using the Dickson CRM method [DAA03]."""
    return olsq(lambda concAcid: solve.DAA03(massAcid, emf, tempK, massSample,
        concAcid, concTotals, eqConstants)['x'][0] - alkCert, 0.1, method='lm')

def halfGran(massAcid, emf, tempK, massSample, alkCert, concTotals, eqConstants):
    """Calibrate acid concentration using half-Gran method [H15]."""
    return olsq(lambda concAcid: solve.halfGran(massAcid, emf, tempK,
        massSample, concAcid, concTotals, eqConstants,
        suppressWarnings=True)[0] - alkCert, 0.1, method='lm')
