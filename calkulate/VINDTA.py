# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019  Matthew Paul Humphreys  (GNU GPLv3)

from . import calib, conc, density, dissoc, io, sim, solve
from numpy import logical_and
from numpy import max as np_max

# ================================================ INPUTS AND THEIR UNITS =====


# Vsamp    = sample volume                    / ml
# Cacid    = acid molality                    / mol/kg
# S        = practical salinity
# AT_cert  = certified total alkalinity       / mol/kg-sw
# CT       = dissolved inorganic carbon       / mol/kg-sw
# PT       = phosphate                        / mol/kg-sw
# SiT      = silicate                         / mol/kg-sw
# Tk_force = titration temperature (optional) / K


# ================================================== PREPARATORY FUNCTION =====


def prep(datfile, Vsamp, S, CT, PT, SiT, burette_cx=1, Tk_force=None):

    Vacid, EMF, Tk = io.VINDTA(datfile)

    if Tk_force is not None:
        Tk[:] = Tk_force

    Msamp = Vsamp * density.sw(Tk[0], S) * 1e-3 # sample mass / kg
    Macid = burette_cx * Vacid * density.acid(Tk) * 1e-3 # acid mass / kg

    XT = conc.XT(S, CT, PT, SiT) # total concentrations  / mol/kg-sw
    KX = dissoc.KX_F(Tk, S, XT[3], XT[4]) # dissociation constants

    return Macid, EMF, Tk, Msamp, XT, KX


# =================================================== HALF-GRAN FUNCTIONS =====


def halfGran(datfile, Vsamp, Cacid, S, CT, PT, burette_cx=1, Tk_force=None):

    Macid, EMF, Tk, Msamp, XT, KX = \
        prep(datfile, Vsamp, S, CT, PT, 0, burette_cx, Tk_force)

    return solve.halfGran(Macid, EMF, Tk, Msamp, Cacid, *XT, *KX)


def halfGran_CRM(datfile, Vsamp, AT_cert, S, CT, PT, burette_cx=1,
    Tk_force=None):

    Macid, EMF, Tk, Msamp, XT, KX = \
        prep(datfile, Vsamp, S, CT, PT, 0, burette_cx, Tk_force)

    Cacid = calib.halfGran(Macid, EMF, Tk, Msamp, AT_cert, XT, KX)['x'][0]

    AT, EMF0, _, _, _, _ = solve.halfGran(Macid, EMF, Tk, Msamp, Cacid,
        *XT, *KX)

    return Cacid, AT, EMF0


# ================================================ PLOT THE LOT FUNCTIONS =====


def guessGran(datfile, Vsamp, Cacid, S, burette_cx=1, Tk_force=None):

    Macid, EMF, Tk, Msamp, _, _ = \
        prep(datfile, Vsamp, S, 0, 0, 0, burette_cx, Tk_force)

    # Evaluate f1 function and corresponding logical
    f1g = solve.f1(Macid, EMF, Tk, Msamp)

    Lg = logical_and(f1g > 0.1 * np_max(f1g), f1g < 0.9 * np_max(f1g))

    # Get first guesses
    ATg, EMF0g, _, pHg = solve.guessGran(Macid, EMF, Tk, Msamp, Cacid)
    EMF0gvec = solve.Gran_EMF0(Macid, EMF, Tk, Msamp, Cacid, ATg)

    # Select data for fitting
    L = logical_and(pHg > 3, pHg < 4)

    return Macid, EMF, Tk, Msamp,  f1g, Lg, EMF0gvec, ATg, EMF0g, pHg, L


def simH(Macid, Tk, Msamp, Cacid, S, AT, CT=0, PT=0, SiT=0):

    XT = conc.XT(S, CT, PT, SiT) # total concentrations  / mol/kg-sw
    XT[0] = AT
    KX = dissoc.KX_F(Tk, S, XT[3], XT[4]) # dissociation constants

    return sim.H(Macid, Msamp, Cacid, XT, KX)


def simAT(Macid,Tk,H,Msamp,S,CT=0,PT=0,SiT=0):

    mu = Msamp / (Msamp + Macid)

    XT = conc.XT(S, CT, PT, SiT) # total concentrations  / mol/kg-sw
    KX = dissoc.KX_F(Tk, S, XT[3], XT[4]) # dissociation constants

    return sim.AT(H, mu, *XT, *KX)


# ========================================================= MPH FUNCTIONS =====


def MPH(datfile, Vsamp, Cacid, S, CT, PT, SiT, burette_cx=1, Tk_force=None):

    Macid, EMF, Tk, Msamp, XT, KX = \
        prep(datfile, Vsamp,S, CT, PT, SiT, burette_cx, Tk_force)

    return solve.MPH(Macid, EMF, Tk, Msamp, Cacid,*XT, KX)


def MPH_CRM(datfile, Vsamp, AT_cert, S, CT, PT, SiT,
    burette_cx=1, Tk_force=None):

    Macid, EMF, Tk, Msamp, XT, KX = \
        prep(datfile, Vsamp,S, CT, PT, SiT, burette_cx, Tk_force)

    Cacid = calib.MPH(Macid, EMF, Tk, Msamp, AT_cert, XT, KX)['x'][0]

    AT, EMF0 = solve.MPH(Macid, EMF, Tk, Msamp, Cacid, *XT, KX)['x']

    return Cacid, AT, EMF0
