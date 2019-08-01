# Calkulate: seawater total alkalinity from titration data.
# Copyright (C) 2019  Matthew Paul Humphreys  (GNU GPLv3)
"""Simulate total alkalinity and pH during titrations."""
from scipy.optimize import least_squares as olsq
from numpy import full_like, nan

def alk(h, mu, XT, KXF, dicFraction=1):
    """Simulate total alkalinity from known pH and total concentrations."""
    OH = KXF['w']/h
    if 'C' in XT.keys():
        co2aq = dicFraction*mu*XT['C']/(1 + KXF['C1']/h +
            KXF['C1']*KXF['C2']/h**2)
        bicarb = KXF['C1']*co2aq/h
        carb = KXF['C2']*bicarb/h
    else:
        co2aq = bicarb = carb = 0
    if 'B' in XT.keys():
        B4 = mu*XT['B']*KXF['B']/(h + KXF['B'])
    else:
        B4 = 0
    if 'S' in XT.keys():
        HSO4 = mu*XT['S']*h/(KXF['S'] + h)
    else:
        HSO4 = 0
    if 'F' in XT.keys():
        HF = mu*XT['F']*h/(KXF['F'] + h)
    else:
        HF = 0
    if 'P' in XT.keys():
        P0 = mu*XT['P']/(1 + KXF['P1']/h + KXF['P1']*KXF['P2']/h**2 +
            KXF['P1']*KXF['P2']*KXF['P3']/h**3)
        P2 = mu*XT['P']/(h**2/(KXF['P1']*KXF['P2']) + h/KXF['P2'] + 1 +
            KXF['P3']/h)
        P3 = mu*XT['P']/(h**3/(KXF['P1']*KXF['P2']*KXF['P3']) +
            h**2/(KXF['P2']*KXF['P3']) + h/KXF['P3'] + 1)
    else:
        P0 = P2 = P3 = 0
    if 'Si' in XT.keys():
        SiOOH3 = mu*XT['Si']*KXF['Si']/(h + KXF['Si'])
    else:
        SiOOH3 = 0
    alk = bicarb + 2*carb + B4 + OH - h - HSO4 - HF - P0 + P2 + 2*P3 + SiOOH3
    return alk, [bicarb, 2*carb, B4, OH, -h, -HSO4, -HF, -P0, P2, 2*P3, SiOOH3]

def pH(massAcid, massSample, concAcid, alk0, XT, KXF):
    """Simulate pH from known total alkalinity and total concentrations."""
    pH = full_like(massAcid, nan)
    mu = massSample/(massSample + massAcid)
    for i, massAcid_i in enumerate(massAcid):
        pH[i] = olsq(lambda pH: alk(10.0**-pH, mu[i], XT, KXF)[0] - mu[i]*alk0 +
                massAcid_i*concAcid/(massAcid_i + massSample),
            8.0, method='lm')['x']
    return pH
