# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019  Matthew Paul Humphreys  (GNU GPLv3)

from numpy import arange, array, full_like, genfromtxt
from .constants import Tzero


def vindta(datfile):

    tdata = genfromtxt(datfile, delimiter='\t', skip_header=2)

    Vacid = tdata[:, 0] # ml
    EMF = tdata[:, 1] # mV
    tempK = tdata[:, 2] + Tzero # K

    return Vacid, EMF, tempK


def Dickson1981():
# Import simulated titration from Dickson (1981) Table 1
#  doi:10.1016/0198-0149(81)90121-7

    # Acid mass in kg
    Macid = arange(0, 2.51, 0.05) * 1e-3

    # Free scale pH
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
    ])

    # Define other variables
    tempK = full_like(Macid, 298.15) # K
    Cacid = 0.3 # mol/kg-soln
    Msamp = 0.2 # kg
    psal = 35.0  # practical

    # Set concentrations
    AT  = 0.00245 # Alkalinity / mol/kg-sw
    BT  = 0.00042 # Borate     / mol/kg-sw
    CT  = 0.00220 # Carbon     / mol/kg-sw
    ST  = 0.02824 # Sulfate    / mol/kg-sw
    FT  = 0.00007 # Fluoride   / mol/kg-sw
    PT  = 0       # Phosphate  / mol/kg-sw
    SiT = 0       # Silicate   / mol/kg-sw
    XT  = [AT, CT, BT, ST, FT, PT, SiT]

    # Set dissociation constants (all are Free pH scale)
    KH2O   = full_like(Macid, 4.32e-14)
    KC1    = full_like(Macid, 1.00e-06)
    KC1KC2 = full_like(Macid, 8.20e-16)
    KC2    = KC1KC2 / KC1
    KB     = full_like(Macid, 1.78e-09)
    KHSO4  = 1 / full_like(Macid, 1.23e+01)
    KHF    = 1 / full_like(Macid, 4.08e+02)
    KP1    = full_like(Macid, 56.8)
    KP2    = full_like(Macid, 8e-7)
    KP3    = full_like(Macid, 1.32e-15) / KP2
    KSi    = full_like(Macid, 1)
    KX = [KC1, KC2, KB, KH2O, KHSO4, KHF, KP1, KP2, KP3, KSi]

    return Macid, pH, tempK, Msamp, Cacid, psal, XT, KX
