# Calkulate: seawater total alkalinity from titration data.
# Copyright (C) 2019  Matthew Paul Humphreys  (GNU GPLv3)
"""Convenient function wrappers for working with VINDTA titration data."""
from . import (calibrate, concentrations, density, dissociation, io, simulate,
    solve)
from numpy import logical_and
from numpy import max as np_max
# ================================================ INPUTS AND THEIR UNITS =====
# volSample = sample volume in ml
# concAcid = acid molality in mol/kg
# pSal = practical salinity (dimensionless)
# alkCert = certified total alkalinity in mol/kg-sw
# CT = dissolved inorganic carbon in mol/kg-sw
# PT = phosphate in mol/kg-sw
# SiT = silicate in mol/kg-sw
# tempKForce = titration temperature (optional) in K

def prep(datfile, volSample, pSal, CT, PT, SiT, buretteCorrection=1,
        tempKForce=None):
    """Import VINDTA-style .dat file and prepare data for analysis."""
    volAcid, emf, tempK = io.vindta(datfile)
    if tempKForce is not None:
        tempK[:] = tempKForce
    massSample = volSample*density.sw(tempK[0], pSal)*1e-3
    massAcid = buretteCorrection*volAcid*density.acid(tempK)*1e-3
    XT = concentrations.XT(pSal, CT, PT, SiT)
    KXF = dissociation.KXF(tempK, pSal, XT)
    return massAcid, emf, tempK, massSample, XT, KXF

# =================================================== HALF-GRAN FUNCTIONS =====
def halfGran(datfile, volSample, concAcid, pSal, CT, PT, buretteCorrection=1,
        tempKForce=None):
    """Solve for alkalinity using the half-Gran method (Humphreys, 2015) for a
    VINDTA-style titration data file.
    """
    massAcid, emf, tempK, massSample, XT, KX = prep(datfile, volSample, pSal,
        CT, PT, 0, buretteCorrection, tempKForce)
    return solve.halfGran(massAcid, emf, tempK, massSample, concAcid, *XT, *KX)

def halfGranCRM(datfile, volSample, alkCert, pSal, CT, PT, buretteCorrection=1,
        tempKForce=None):
    """Solve for acid concentration using the half-Gran method (Humphreys, 2015)
    for a VINDTA-style titration data file.
    """
    massAcid, emf, tempK, massSample, XT, KX = prep(datfile, volSample, pSal,
        CT, PT, 0, buretteCorrection, tempKForce)
    concAcid = calibrate.halfGran(massAcid, emf, tempK, massSample, alkCert,
        XT, KX)['x'][0]
    alk, emf0, _, _, _, _ = solve.halfGran(massAcid, emf, tempK, massSample,
        concAcid, *XT, *KX)
    return concAcid, alk, emf0

# ================================================ PLOT THE LOT FUNCTIONS =====
def guessGran(datfile, volSample, concAcid, pSal, buretteCorrection=1,
        tempKForce=None):
    massAcid, EMF, tempK, massSample, _, _ = prep(datfile, volSample, pSal, 0,
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

def simH(massAcid, tempK, massSample, concAcid, pSal, alk, CT=0, PT=0, SiT=0):
    XT = concentrations.XT(pSal, CT, PT, SiT)
    XT[0] = alk
    KX = dissociation.KXF(tempK, pSal, XT)
    return simulate.H(massAcid, massSample, concAcid, XT, KX)

def simAT(massAcid, tempK, H, massSample, pSal, CT=0, PT=0, SiT=0):
    mu = solve.mu(massAcid, massSample)
    XT = concentrations.XT(pSal, CT, PT, SiT)
    KX = dissociation.KXF(tempK, pSal, XT)
    return simulate.alk(H, mu, XT, KX)

def complete(datfile, volSample, concAcid, pSal, CT, PT, SiT,
        buretteCorrection=1, tempKForce=None):
    """Solve for alkalinity using the Complete Calculation method for a
    VINDTA-style titration data file.
    """
    massAcid, emf, tempK, massSample, XT, KX = prep(datfile, volSample, pSal,
        CT, PT, SiT, buretteCorrection, tempKForce)
    return solve.complete(massAcid, emf, tempK, massSample, concAcid, XT, KX)

def completeCRM(datfile, volSample, alkCert, pSal, CT, PT, SiT,
        buretteCorrection=1, tempKForce=None):
    """Solve for acid concentration using the Complete Calculation method for a
    VINDTA-style titration data file.
    """
    massAcid, emf, tempK, massSample, XT, KX = prep(datfile, volSample, pSal,
        CT, PT, SiT, buretteCorrection, tempKForce)
    concAcid = calibrate.complete(massAcid, emf, tempK, massSample, alkCert,
        XT, KX)['x'][0]
    alk, emf0 = solve.complete(massAcid, emf, tempK, massSample, concAcid, XT,
        KX)['x']
    return concAcid, alk, emf0
