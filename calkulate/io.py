# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019  Matthew Paul Humphreys  (GNU GPLv3)

from numpy import arange, array, full_like, genfromtxt
from .constants import Tzero

def vindta(datfile):
    """Import VINDTA-style .dat file titration table."""
    tdata = genfromtxt(datfile, delimiter='\t', skip_header=2)
    Vacid = tdata[:, 0] # ml
    EMF = tdata[:, 1] # mV
    tempK = tdata[:, 2] + Tzero # K
    return Vacid, EMF, tempK

def Dickson1981():
    """Import simulated titration from Dickson (1981), Table 1."""
    Macid = arange(0, 2.51, 0.05) * 1e-3 # Acid mass in kg
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
    tempK = full_like(Macid, 298.15) # K
    Cacid = 0.3 # mol/kg-soln
    Msamp = 0.2 # kg
    psal = 35.0 # practical salinity
    # Set concentrations, all in mol/kg-sw
    AT = 0.00245 # Alkalinity
    XT = {}
    XT['B'] = 0.00042 # Borate
    XT['C'] = 0.00220 # Carbon
    XT['S'] = 0.02824 # Sulfate
    XT['F'] = 0.00007 # Fluoride
    XT['P'] = 0 # Phosphate
    XT['Si'] = 0 # Silicate
    # Set dissociation constants, all on Free pH scale
    KXF = {}
    KXF['w'] = full_like(Macid, 4.32e-14)
    KXF['C1'] = full_like(Macid, 1.00e-06)
    KC1KC2 = full_like(Macid, 8.20e-16)
    KXF['C2'] = KC1KC2 / KXF['C1']
    KXF['B'] = full_like(Macid, 1.78e-09)
    KXF['S'] = 1 / full_like(Macid, 1.23e+01)
    KXF['F'] = 1 / full_like(Macid, 4.08e+02)
    KXF['P1'] = full_like(Macid, 56.8)
    KXF['P2'] = full_like(Macid, 8e-7)
    KXF['P3'] = full_like(Macid, 1.32e-15) / KXF['P2']
    KXF['Si'] = full_like(Macid, 1)
    return Macid, pH, tempK, Msamp, Cacid, psal, AT, XT, KXF
