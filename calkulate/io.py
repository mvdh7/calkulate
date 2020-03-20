# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew Paul Humphreys  (GNU GPLv3)
"""Import data to use with Calkulate and export results."""
from numpy import arange, array, full_like, genfromtxt
from .constants import tZero

def datfile(datFile, delimiter='\t', skip_header=2, **kwargs):
    """Import a single VINDTA-style .dat file titration dataset."""
    tData = genfromtxt(datFile, delimiter=delimiter, skip_header=skip_header,
        **kwargs)
    volAcid = tData[:, 0] # ml
    emf = tData[:, 1] # mV
    tempK = tData[:, 2] + tZero # K
    return volAcid, emf, tempK

# Add alias to avoid breaking old code
vindta = datfile

def writeDat(datFile, volAcid, emf, tempK, line0='', line1=''):
    """Write a titration dataset to a VINDTA-style .dat file."""
    with open(datFile, 'w') as f:
        f.write('{}\n{}\n'.format(line0, line1))
        for i in range(len(volAcid)):
            f.write('{:.5f}\t{:.5f}\t{:.3f}\n'.format(
                volAcid[i], emf[i], tempK[i]-tZero))

def Dickson1981(withPhosphate=True):
    """Import simulated titrations from D81."""
    massAcid = arange(0, 2.51, 0.05)*1e-3 # acid mass in kg
    tempK = full_like(massAcid, 298.15) # K
    concAcid = 0.3 # mol/kg-soln
    massSample = 0.2 # kg
    pSal = 35.0 # practical salinity
    # Set concentrations, all in mol/kg-sw
    alk = 0.00245 # Alkalinity
    concTotals = {}
    concTotals['B'] = 0.00042 # Borate
    concTotals['C'] = 0.00220 # Carbon
    concTotals['S'] = 0.02824 # Sulfate
    concTotals['F'] = 0.00007 # Fluoride
    concTotals['Si'] = 0 # Silicate
    # Set dissociation constants, all on Free pH scale
    eqConstants = {}
    eqConstants['w'] = full_like(massAcid, 4.32e-14)
    eqConstants['C1'] = full_like(massAcid, 1.00e-06)
    kC1kC2 = full_like(massAcid, 8.20e-16)
    eqConstants['C2'] = kC1kC2/eqConstants['C1']
    eqConstants['B'] = full_like(massAcid, 1.78e-09)
    eqConstants['S'] = 1/full_like(massAcid, 1.23e+01)
    eqConstants['F'] = 1/full_like(massAcid, 4.08e+02)
    eqConstants['P1'] = full_like(massAcid, 1/56.8)
    eqConstants['P2'] = full_like(massAcid, 8e-7)
    eqConstants['P3'] = full_like(massAcid, 1.32e-15)/eqConstants['P2']
    eqConstants['Si'] = full_like(massAcid, 1)
    if withPhosphate:
        concTotals['P'] = 0.00001 # Phosphate
        pH = array([
            8.046129, 7.902742, 7.724342, 7.509157, 7.284827, 7.089577,
            6.931664, 6.802646, 6.693843, 6.588221, 6.514787, 6.437841,
            6.366846, 6.299338, 6.235339, 6.173652, 6.113581, 6.054522,
            5.995928, 5.937273, 5.878028, 5.817630, 5.755453, 5.690759,
            5.622646, 5.549957, 5.471112, 5.383919, 5.285082, 5.169392,
            5.028069, 4.845783, 4.600713, 4.304478, 4.046487, 3.860077,
            3.723045, 3.616672, 3.530308, 3.457817, 3.395443, 3.340748,
            3.292070, 3.248228, 3.208353, 3.171793, 3.138041, 3.106699,
            3.077446, 3.050022, 3.024213,
        ]) # Free scale pH
    else:
        concTotals['P'] = 0 # Phosphate
        pH = array([
            8.065650, 7.925895, 7.752062, 7.539922, 7.312923, 7.111723,
            6.948715, 6.816116, 6.704825, 6.608421, 6.522665, 6.444704,
            6.372551, 6.304758, 6.240231, 6.178101, 6.117655, 6.058275,
            5.999403, 5.940504, 5.881044, 5.820455, 5.758107, 5.693259,
            5.625006, 5.552183, 5.473224, 5.385920, 5.286975, 5.171175,
            5.029724, 4.847252, 4.601818, 4.304988, 4.046597, 3.860034,
            3.722943, 3.616542, 3.530163, 3.457664, 3.395285, 3.340587,
            3.291906, 3.248062, 3.208187, 3.171626, 3.137874, 3.106531,
            3.077278, 3.049854, 3.024045,
        ]) # Free scale pH
    return (massAcid, pH, tempK, massSample, concAcid, pSal, alk, concTotals,
        eqConstants)
