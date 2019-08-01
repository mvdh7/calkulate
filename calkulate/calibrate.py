# Calkulate: seawater total alkalinity from titration data.
# Copyright (C) 2019  Matthew Paul Humphreys  (GNU GPLv3)
"""Calibrate acid titrant concentration using certified reference material
measurements.
"""
from scipy.optimize import least_squares as olsq
from . import solve

def halfGran(massAcid, emf, tempK, massSample, alkCert, XT, KXF):
    """Calibrate acid concentration using half-Gran method."""
    return olsq(lambda concAcid: solve.halfGran(massAcid, emf, tempK,
        massSample, concAcid, XT, KXF, suppressWarnings=True)[0] - alkCert,
        0.1, method='lm')

def complete(massAcid, emf, tempK, massSample, alkCert, XT, KXF):
    """Calibrate acid concentration using the Complete Calculation."""
    return olsq(lambda concAcid: solve.complete(massAcid, emf, tempK,
        massSample, concAcid, XT, KXF)['x'][0] - alkCert, 0.1, method='lm')

def DAA03(massAcid, emf, tempK, massSample, alkCert, XT, KXF):
    """Calibrate acid concentration following Dickson et al. (2003)."""
    return olsq(lambda concAcid: solve.DAA03(massAcid, emf, tempK, massSample,
        concAcid, XT, KXF)['x'][0] - alkCert, 0.1, method='lm')
