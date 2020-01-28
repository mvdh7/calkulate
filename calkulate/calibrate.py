# Calkulate: seawater total alkalinity from titration data.
# Copyright (C) 2019-2020  Matthew Paul Humphreys  (GNU GPLv3)
"""Calibrate acid titrant concentration using certified reference material
measurements.
"""
from scipy.optimize import least_squares as olsq
from . import solve

def concAcid(massAcid, emf, tempK, massSample, alkCert, concTotals,
        eqConstants, solver='complete', **kwargs):
    """Calibrate the acid concentration using known sample alkalinity."""
    if solver.lower() in solve.allSolvers.keys():
        solveFunc = solve.allSolvers[solver.lower()]
        concAcidOptResult = olsq(lambda concAcid: solveFunc(massAcid, emf,
            tempK, massSample, concAcid, concTotals, eqConstants,
            **kwargs)['x'][0] - alkCert, 0.1, method='lm')
    else:
        print('calkulate.calibrate.concAcid: solver not recognised.')
        print('Options (case-insensitive):' +
            (len(solve.allSolvers.keys())*' \'{}\'').format(
                *solve.allSolvers.keys()))
        concAcidOptResult = {'x': [None,]}
    return concAcidOptResult
