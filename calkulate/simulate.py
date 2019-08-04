# Calkulate: seawater total alkalinity from titration data.
# Copyright (C) 2019  Matthew Paul Humphreys  (GNU GPLv3)
"""Simulate total alkalinity and pH during titrations."""
from scipy.optimize import least_squares as olsq
from numpy import full_like, nan
from . import solve

def alk(h, mu, concTotals, eqConstants):
    """Simulate total alkalinity from known pH and total concentrations."""
    OH = eqConstants['w']/h
    if 'C' in concTotals.keys():
        co2aq = mu*concTotals['C']/(1 + eqConstants['C1']/h +
            eqConstants['C1']*eqConstants['C2']/h**2)
        bicarb = eqConstants['C1']*co2aq/h
        carb = eqConstants['C2']*bicarb/h
    else:
        bicarb = carb = 0
    if 'B' in concTotals.keys():
        B4 = mu*concTotals['B']*eqConstants['B']/(h + eqConstants['B'])
    else:
        B4 = 0
    if 'S' in concTotals.keys():
        HSO4 = mu*concTotals['S']*h/(eqConstants['S'] + h)
    else:
        HSO4 = 0
    if 'F' in concTotals.keys():
        HF = mu*concTotals['F']*h/(eqConstants['F'] + h)
    else:
        HF = 0
    if 'P' in concTotals.keys():
        P0 = mu*concTotals['P']/(1 + eqConstants['P1']/h +
            eqConstants['P1']*eqConstants['P2']/h**2 +
            eqConstants['P1']*eqConstants['P2']*eqConstants['P3']/h**3)
        P2 = mu*concTotals['P']/(h**2/(eqConstants['P1']*eqConstants['P2']) +
            h/eqConstants['P2'] + 1 + eqConstants['P3']/h)
        P3 = mu*concTotals['P']/(h**3/(eqConstants['P1']*eqConstants['P2']*
            eqConstants['P3']) + h**2/(eqConstants['P2']*eqConstants['P3']) +
            h/eqConstants['P3'] + 1)
    else:
        P0 = P2 = P3 = 0
    if 'Si' in concTotals.keys():
        SiOOH3 = mu*concTotals['Si']*eqConstants['Si']/(h + eqConstants['Si'])
    else:
        SiOOH3 = 0
    alk = bicarb + 2*carb + B4 + OH - h - HSO4 - HF - P0 + P2 + 2*P3 + SiOOH3
    components = {
        '+HCO3': bicarb,
        '+2*CO3': 2*carb,
        '+B(OH)4': B4,
        '+OH': OH,
        '-H': -h,
        '-HSO4': -HSO4,
        '-HF': -HF,
        '-H3PO4': -P0,
        '+HPO4': P2,
        '+2*PO4': P3,
        '+SiO(OH)3': SiOOH3,
    }
    return alk, components

def pH(massAcid, massSample, concAcid, alk0, concTotals, eqConstants):
    """Simulate pH from known total alkalinity and total concentrations."""
    pH = full_like(massAcid, nan)
    mu = solve.mu(massAcid, massSample)
    for i, massAcid_i in enumerate(massAcid):
        pH[i] = olsq(lambda pH: alk(10.0**-pH, mu[i], concTotals,
                eqConstants)[0] - mu[i]*alk0 + massAcid_i*concAcid/
                (massAcid_i + massSample),
            8.0, method='lm')['x'][0]
    return pH
