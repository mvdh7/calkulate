# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019  Matthew Paul Humphreys  (GNU GPLv3)

from scipy.optimize import least_squares as olsq
from numpy import full_like, inf, nan

def AT(H, mu, XT, KXF):
    """Simulate total alkalinity from known pH and total concentrations."""
    CO2aq = mu * XT['C'] / (1 + KXF['C1']/H + KXF['C1']*KXF['C2']/H**2)
    bicarb = KXF['C1'] * CO2aq / H
    carb = KXF['C2'] * bicarb / H
    B4 = mu * XT['B']*KXF['B'] / (H + KXF['B'])
    OH = KXF['w'] / H
    HSO4 = mu * XT['S']*H / (KXF['S'] + H)
    HF = mu * XT['F']*H / (KXF['F'] + H)
    P0 = mu * XT['P'] / (1 + KXF['P1']/H + KXF['P1']*KXF['P2']/H**2 \
        + KXF['P1']*KXF['P2']*KXF['P3']/H**3)
    P2 = mu * XT['P'] / (H**2/(KXF['P1']*KXF['P2']) + H/KXF['P2'] + 1 \
        + KXF['P3']/H)
    P3 = mu * XT['P'] / (H**3/(KXF['P1']*KXF['P2']*KXF['P3']) \
        + H**2/(KXF['P2']*KXF['P3']) + H/KXF['P3'] + 1)
    SiOOH3 = mu * XT['Si']*KXF['Si'] / (H + KXF['Si'])
    AT = bicarb + 2*carb + B4 + OH - H - HSO4 - HF - P0 + P2 + 2*P3 + SiOOH3
    return AT, [bicarb, 2*carb, B4, OH, -H, -HSO4, -HF, -P0, P2, 2*P3, SiOOH3]

def H(Macid, Msamp, Cacid, AT, XT, KXF):
    """Simulate pH from known total alkalinity and total concentrations."""
    H = full_like(Macid, nan)
    mu = Msamp / (Msamp + Macid)
    for i, Macid_i in enumerate(Macid):
        H[i] = olsq(lambda H: \
            AT(H, mu[i], XT, KXF)[0] - mu[i]*AT + \
                Macid_i*Cacid/(Macid_i + Msamp),
            1e-8, method='trf', bounds=(1e-15, inf))['x']
    return H
