# Calkulate: seawater total alkalinity from titration data.
# Copyright (C) 2019  Matthew Paul Humphreys  (GNU GPLv3)
"""Calibrate acid titrant concentration using certified reference material
measurements.
"""
from scipy.optimize import least_squares as olsq
from . import solve

def concAcid(massAcid, emf, tempK, massSample, alkCert, concTotals,
        eqConstants, solver='complete', **kwargs):
    """Calibrate the acid concentration using known sample alkalinity."""
    solveFuncs = {
        'complete': solve.complete,
        'daa03': solve.DAA03,
        'dickson1981': solve.Dickson1981,
        'halfgran': solve.halfGran,
    }
    if solver.lower() in solveFuncs.keys():
        solveFunc = solveFuncs[solver.lower()]
        concAcidOptResult = olsq(lambda concAcid: solveFunc(massAcid, emf,
            tempK, massSample, concAcid, concTotals, eqConstants,
            **kwargs)['x'][0] - alkCert, 0.1, method='lm')
    else:
        print('calkulate.calibrate: solver not recognised.')
        print('Options (case-insensitive):' +
            (len(solveFuncs.keys())*' \'{}\'').format(*solveFuncs.keys()))
        concAcidOptResult = {'x': [None,]}
    return concAcidOptResult
