# Calkulate: seawater total alkalinity from titration data.
# Copyright (C) 2019-2020  Matthew Paul Humphreys  (GNU GPLv3)
from . import density

def vol2massAcid(volAcid, tempK, buretteCorrection=1.0):
    """Convert acid volume(s) in ml to mass(es) in kg."""
    return buretteCorrection*volAcid*density.acid(tempK)*1e-3

def vol2massSample(volSample, tempK, pSal):
    """Convert sample volume in ml to mass in kg."""
    return volSample*density.sw(tempK, pSal)*1e-3
