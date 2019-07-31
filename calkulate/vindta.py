# Calkulate: seawater total alkalinity from titration data.
# Copyright (C) 2019  Matthew Paul Humphreys  (GNU GPLv3)
"""Convenient function wrappers for working with VINDTA titration data."""
from . import (calibrate, concentrations, density, dissociation, io,
    simulate, solve)
from numpy import logical_and
from numpy import max as np_max

# ================================================ INPUTS AND THEIR UNITS =====

# Vsamp      = sample volume                    in ml
# Cacid      = acid molality                    in mol/kg
# S          = practical salinity
# AT_cert    = certified total alkalinity       in mol/kg-sw
# CT         = dissolved inorganic carbon       in mol/kg-sw
# PT         = phosphate                        in mol/kg-sw
# SiT        = silicate                         in mol/kg-sw
# tempKforce = titration temperature (optional) in K

# ================================================== PREPARATORY FUNCTION =====

def prep(datfile, Vsamp, psal, CT, PT, SiT, burette_cx=1, tempKforce=None):
    """Import VINDTA-style .dat file and prepare data for analysis."""
    Vacid, EMF, tempK = io.vindta(datfile)
    if tempKforce is not None:
        tempK[:] = tempKforce
    Msamp = Vsamp * density.sw(tempK[0], psal) * 1e-3 # sample mass / kg
    Macid = burette_cx * Vacid * density.acid(tempK) * 1e-3 # acid mass / kg
    XT = concentrations.XT(psal, CT, PT, SiT)
    KXF = dissociation.KXF(tempK, psal, XT)
    return Macid, EMF, tempK, Msamp, XT, KXF


# =================================================== HALF-GRAN FUNCTIONS =====

def halfGran(datfile, Vsamp, Cacid, psal, CT, PT, burette_cx=1,
        tempKforce=None):
    Macid, EMF, tempK, Msamp, XT, KX = prep(datfile, Vsamp, psal, CT, PT, 0, 
        burette_cx, tempKforce)
    return solve.halfGran(Macid, EMF, tempK, Msamp, Cacid, *XT, *KX)

def halfGranCRM(datfile, Vsamp, AT_cert, psal, CT, PT, burette_cx=1,
        tempKforce=None):
    Macid, EMF, tempK, Msamp, XT, KX = prep(datfile, Vsamp, psal, CT, PT, 0,
        burette_cx, tempKforce)
    Cacid = calibrate.halfGran(Macid, EMF, tempK, Msamp, AT_cert,
        XT, KX)['x'][0]
    AT, EMF0, _, _, _, _ = solve.halfGran(Macid, EMF, tempK, Msamp, Cacid,
        *XT, *KX)
    return Cacid, AT, EMF0


# ================================================ PLOT THE LOT FUNCTIONS =====

def guessGran(datfile, Vsamp, Cacid, psal, burette_cx=1, tempKforce=None):
    Macid, EMF, tempK, Msamp, _, _ = prep(datfile, Vsamp, psal, 0, 0, 0,
        burette_cx, tempKforce)
    # Evaluate f1 function and corresponding logical
    f1g = solve.f1(Macid, EMF, tempK, Msamp)
    Lg = logical_and(f1g > 0.1 * np_max(f1g), f1g < 0.9 * np_max(f1g))
    # Get first guesses
    ATg, EMF0g, _, pHg = solve.guessGran(Macid, EMF, tempK, Msamp, Cacid)
    EMF0gvec = solve.Gran_EMF0(Macid, EMF, tempK, Msamp, Cacid, ATg)
    # Select data for fitting
    L = logical_and(pHg > 3, pHg < 4)
    return Macid, EMF, tempK, Msamp,  f1g, Lg, EMF0gvec, ATg, EMF0g, pHg, L

def simH(Macid, tempK, Msamp, Cacid, psal, AT, CT=0, PT=0, SiT=0):
    XT = concentrations.XT(psal, CT, PT, SiT) # total concentrations
    XT[0] = AT
    KX = dissociation.KXF(tempK, psal, XT) # dissociation constants
    return simulate.H(Macid, Msamp, Cacid, XT, KX)

def simAT(Macid, tempK, H, Msamp, psal, CT=0, PT=0, SiT=0):
    mu = Msamp / (Msamp + Macid)
    XT = concentrations.XT(psal, CT, PT, SiT) # total concentrations
    KX = dissociation.KXF(tempK, psal, XT) # dissociation constants
    return simulate.AT(H, mu, XT, KX)


# ========================================================= MPH FUNCTIONS =====

def complete(datfile, Vsamp, Cacid, psal, CT, PT, SiT, burette_cx=1, 
        tempKforce=None):
    Macid, EMF, tempK, Msamp, XT, KX = prep(datfile, Vsamp, psal, CT, PT, SiT,
        burette_cx, tempKforce)
    return solve.complete(Macid, EMF, tempK, Msamp, Cacid, XT, KX)

def completeCRM(datfile, Vsamp, AT_cert, psal, CT, PT, SiT,
        burette_cx=1, tempKforce=None):
    Macid, EMF, tempK, Msamp, XT, KX = prep(datfile, Vsamp, psal, CT, PT, SiT,
        burette_cx, tempKforce)
    Cacid = calibrate.complete(Macid, EMF, tempK, Msamp, AT_cert,
        XT, KX)['x'][0]
    AT, EMF0 = solve.complete(Macid, EMF, tempK, Msamp, Cacid, XT, KX)['x']
    return Cacid, AT, EMF0
