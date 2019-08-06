# Calkulate: seawater total alkalinity from titration data.
# Copyright (C) 2019  Matthew Paul Humphreys  (GNU GPLv3)
"""Function wrappers for working with VINDTA-style titration data."""
from . import calibrate, concentrations, density, dissociation, io, solve
# ================================================ INPUTS AND THEIR UNITS =====
# datFile = .dat file name (and path)
# volSample = sample volume in ml
# concAcid = acid molality in mol/kg
# pSal = practical salinity (dimensionless)
# alkCert = certified total alkalinity in mol/kg-sw
# totalCarbonate = dissolved inorganic carbon in mol/kg-sw
# totalPhosphate = phosphate in mol/kg-sw
# totalSilicate = silicate in mol/kg-sw
# tempKForce = titration temperature (optional) in K

def prep(datFile, volSample, pSal, totalCarbonate, totalPhosphate,
        totalSilicate, buretteCorrection=1, tempKForce=None):
    """Import VINDTA-style .dat file and prepare data for analysis."""
    volAcid, emf, tempK = io.vindta(datFile)
    if tempKForce is not None:
        tempK[:] = tempKForce
    massSample = volSample*density.sw(tempK[0], pSal)*1e-3
    massAcid = buretteCorrection*volAcid*density.acid(tempK)*1e-3
    concTotals = concentrations.concTotals(pSal, totalCarbonate,
        totalPhosphate, totalSilicate)
    eqConstants = dissociation.eqConstants(tempK, pSal, concTotals)
    return massAcid, emf, tempK, massSample, concTotals, eqConstants

def alk(datFile, volSample, concAcid, pSal, totalCarbonate, totalPhosphate,
        totalSilicate, solver='complete', buretteCorrection=1,
        tempKForce=None, **kwargs):
    """Solve for alkalinity from a VINDTA-style titration .dat file."""
    massAcid, emf, tempK, massSample, concTotals, eqConstants = prep(datFile,
        volSample, pSal, totalCarbonate, totalPhosphate, totalSilicate,
        buretteCorrection, tempKForce)
    if solver.lower() in solve.allSolvers.keys():
        solveFunc = solve.allSolvers[solver.lower()]
        alkOptResult = solveFunc(massAcid, emf, tempK, massSample, concAcid,
            concTotals, eqConstants, **kwargs)
    else:
        print('calkulate.vindta.alk: solver not recognised.')
        print('Options (case-insensitive):' +
            (len(solve.allSolvers.keys())*' \'{}\'').format(
                *solve.allSolvers.keys()))
        alkOptResult = {'x': [None,]}
    return alkOptResult

def concAcid(datFile, volSample, alkCert, pSal, totalCarbonate,
        totalPhosphate, totalSilicate, solver='complete', buretteCorrection=1,
        tempKForce=None, **kwargs):
    """Solve for acid concentration from a VINDTA-style titration .dat file."""
    massAcid, emf, tempK, massSample, concTotals, eqConstants = prep(datFile,
        volSample, pSal, totalCarbonate, totalPhosphate, totalSilicate,
        buretteCorrection, tempKForce)
    concAcidOptResult = calibrate.concAcid(massAcid, emf, tempK, massSample,
        alkCert, concTotals, eqConstants, solver, **kwargs)
    return concAcidOptResult
