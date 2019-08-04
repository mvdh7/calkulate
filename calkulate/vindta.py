# Calkulate: seawater total alkalinity from titration data.
# Copyright (C) 2019  Matthew Paul Humphreys  (GNU GPLv3)
"""Function wrappers for working with VINDTA-style titration data."""
from . import (calibrate, concentrations, density, dissociation, io, simulate,
    solve)
from numpy import logical_and
from numpy import max as np_max
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

# =================================================== HALF-GRAN FUNCTIONS =====
def halfGran(datFile, volSample, concAcid, pSal, totalCarbonate,
        totalPhosphate, buretteCorrection=1, tempKForce=None):
    """Solve for alkalinity using the half-Gran method [H15] for a VINDTA-style
    titration data file.
    """
    massAcid, emf, tempK, massSample, concTotals, eqConstants = prep(datFile,
        volSample, pSal, totalCarbonate, totalPhosphate, 0,
        buretteCorrection, tempKForce)
    return solve.halfGran(massAcid, emf, tempK, massSample, concAcid,
        concTotals, eqConstants)

def halfGranCRM(datFile, volSample, alkCert, pSal, totalCarbonate,
        totalPhosphate, buretteCorrection=1, tempKForce=None):
    """Solve for acid concentration using the half-Gran method [H15]
    for a VINDTA-style titration data file.
    """
    massAcid, emf, tempK, massSample, concTotals, eqConstants = prep(datFile,
        volSample, pSal, totalCarbonate, totalPhosphate, 0, buretteCorrection,
        tempKForce)
    concAcid = calibrate.halfGran(massAcid, emf, tempK, massSample, alkCert,
        concTotals, eqConstants)['x'][0]
    alk, emf0, _, _, _, _ = solve.halfGran(massAcid, emf, tempK, massSample,
        concAcid, *concTotals, *eqConstants)
    return concAcid, alk, emf0

# ================================================ PLOT THE LOT FUNCTIONS =====
def guessGran(datFile, volSample, concAcid, pSal, buretteCorrection=1,
        tempKForce=None):
    massAcid, EMF, tempK, massSample, _, _ = prep(datFile, volSample, pSal, 0,
        0, 0, buretteCorrection, tempKForce)
    # Evaluate f1 function and corresponding logical
    f1g = solve.f1(massAcid, EMF, tempK, massSample)
    Lg = logical_and(f1g > 0.1 * np_max(f1g), f1g < 0.9 * np_max(f1g))
    # Get first guesses
    ATg, EMF0g, _, pHg = solve.guessGran(massAcid, EMF, tempK, massSample,
        concAcid)
    EMF0gvec = solve.Gran_EMF0(massAcid, EMF, tempK, massSample, concAcid, ATg)
    # Select data for fitting
    L = logical_and(pHg > 3, pHg < 4)
    return (massAcid, EMF, tempK, massSample,  f1g, Lg, EMF0gvec, ATg, EMF0g,
        pHg, L)

def simH(massAcid, tempK, massSample, concAcid, pSal, alk, totalCarbonate=0,
        totalPhosphate=0, totalSilicate=0):
    concTotals = concentrations.concTotals(pSal, totalCarbonate,
        totalPhosphate, totalSilicate)
    concTotals[0] = alk
    eqConstants = dissociation.eqConstants(tempK, pSal, concTotals)
    return simulate.H(massAcid, massSample, concAcid, concTotals, eqConstants)

def simAT(massAcid, tempK, H, massSample, pSal, totalCarbonate=0,
        totalPhosphate=0, totalSilicate=0):
    mu = solve.mu(massAcid, massSample)
    concTotals = concentrations.concTotals(pSal, totalCarbonate,
        totalPhosphate, totalSilicate)
    eqConstants = dissociation.eqConstants(tempK, pSal, concTotals)
    return simulate.alk(H, mu, concTotals, eqConstants)

def complete(datFile, volSample, concAcid, pSal, totalCarbonate,
        totalPhosphate, totalSilicate, buretteCorrection=1, tempKForce=None):
    """Solve for alkalinity using the complete calculation method for a
    VINDTA-style titration data file.
    """
    massAcid, emf, tempK, massSample, concTotals, eqConstants = prep(datFile,
        volSample, pSal, totalCarbonate, totalPhosphate, totalSilicate,
        buretteCorrection, tempKForce)
    return solve.complete(massAcid, emf, tempK, massSample, concAcid,
        concTotals, eqConstants)

def completeCRM(datFile, volSample, alkCert, pSal, totalCarbonate,
        totalPhosphate, totalSilicate, buretteCorrection=1, tempKForce=None):
    """Solve for acid concentration using the complete calculation method for a
    VINDTA-style titration data file.
    """
    massAcid, emf, tempK, massSample, concTotals, eqConstants = prep(datFile,
        volSample, pSal, totalCarbonate, totalPhosphate, totalSilicate,
        buretteCorrection, tempKForce)
    concAcid = calibrate.complete(massAcid, emf, tempK, massSample, alkCert,
        concTotals, eqConstants)['x'][0]
    alk, emf0 = solve.complete(massAcid, emf, tempK, massSample, concAcid,
        concTotals, eqConstants)['x']
    return concAcid, alk, emf0
