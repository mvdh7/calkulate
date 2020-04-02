# Calkulate: seawater total alkalinity from titration data.
# Copyright (C) 2019-2020  Matthew Paul Humphreys  (GNU GPLv3)
from numpy import exp, log
from . import density
from .constants import F, R

def vol2massAcid(volAcid, tempK, buretteCorrection=1.0):
    """Convert acid volume(s) in ml to mass(es) in kg."""
    return buretteCorrection*volAcid*density.acid(tempK)*1e-3

def vol2massSample(volSample, tempK, pSal):
    """Convert sample volume in ml to mass in kg."""
    return volSample*density.sw(tempK, pSal)*1e-3

def emf2h(emf, emf0, tempK):
    """Convert EMF to [H+]."""
    # DAA03 Eq. (13) with typo corrected (i.e. EMF and EMF0 switched)
    return exp((emf - emf0)*F/(R*tempK))

def h2emf(h, emf0, tempK):
    """Convert [H+] to EMF."""
    return emf0 + log(h)*R*tempK/F

def f2dEmf0(tempK, f):
    return log(f)*R*tempK/F
