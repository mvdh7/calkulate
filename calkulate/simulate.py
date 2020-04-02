# Calkulate: seawater total alkalinity from titration data.
# Copyright (C) 2019-2020  Matthew Paul Humphreys  (GNU GPLv3)
"""Simulate total alkalinity and pH during titrations."""
from scipy.optimize import least_squares as olsq
from numpy import arange, full_like, size, nan
from . import concentrations, density, dissociation, solve
from . import titration as titr # to avoid namespace clash while maintaining
                                # backwards compatibility

def alk(h, mu, concTotals, eqConstants):
    """Simulate total alkalinity from known pH and total concentrations."""
    # Water
    OH = eqConstants['w']/h
    # Carbonate and bicarbonate
    if 'C' in concTotals.keys():
        KC1 = eqConstants['C1']
        KC2 = eqConstants['C2']
        TC = concTotals['C']
        co2aq = mu*TC/(1 + KC1/h + KC1*KC2/h**2)
        bicarb = KC1*co2aq/h
        carb = KC2*bicarb/h
    else:
        bicarb = carb = 0
    # Borate
    if 'B' in concTotals.keys():
        B4 = mu*concTotals['B']*eqConstants['B']/(h + eqConstants['B'])
    else:
        B4 = 0
    # Bisulfate
    if 'S' in concTotals.keys():
        HSO4 = mu*concTotals['S']*h/(eqConstants['S'] + h)
    else:
        HSO4 = 0
    # Hydrogen fluoride
    if 'F' in concTotals.keys():
        HF = mu*concTotals['F']*h/(eqConstants['F'] + h)
    else:
        HF = 0
    # Phosphates
    if 'P' in concTotals.keys():
        TP = concTotals['P']
        KP1 = eqConstants['P1']
        KP2 = eqConstants['P2']
        KP3 = eqConstants['P3']
        P0 = mu*TP/(1 + KP1/h + KP1*KP2/h**2 + KP1*KP1*KP3/h**3)
        P2 = mu*TP/(h**2/(KP1*KP2) + h/KP2 + 1 + KP3/h)
        P3 = mu*TP/(h**3/(KP1*KP2*KP3) + h**2/(KP2*KP3) + h/KP3 + 1)
    else:
        P0 = P2 = P3 = 0
    # Silicate
    if 'Si' in concTotals.keys():
        SiOOH3 = mu*concTotals['Si']*eqConstants['Si']/(h + eqConstants['Si'])
    else:
        SiOOH3 = 0
    # Ammonia
    if 'NH3' in concTotals.keys():
        TNH3 = concTotals['NH3']
        KNH3 = eqConstants['NH3']
        NH3 = TNH3*KNH3/(KNH3 + h)
    else:
        NH3 = 0
    # Hydrogen sulfide
    if 'H2S' in concTotals.keys():
        TH2S = concTotals['H2S']
        KH2S = concTotals['H2S']
        HS = TH2S*KH2S/(KH2S + h)
    else:
        HS = 0
    alk = (bicarb + 2*carb + B4 + OH - h - HSO4 - HF - P0 + P2 + 2*P3 +
        SiOOH3 + NH3 + HS)
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
        '+NH3': NH3,
        '+HS': HS,
    }
    return alk, components

def pH(massAcid, massSample, concAcid, alk0, concTotals, eqConstants):
    """Simulate pH from known total alkalinity and total concentrations."""
    pH = full_like(massAcid, nan)
    mu = solve.mu(massAcid, massSample)
    for i, massAcid_i in enumerate(massAcid):
        iConcTotals = {k: v if size(v)==1 else v[i]
            for k, v in concTotals.items()}
        pH[i] = olsq(lambda pH: alk(10.0**-pH, mu[i], iConcTotals,
                eqConstants)[0] - mu[i]*alk0 + massAcid_i*concAcid/
                (massAcid_i + massSample),
            8.0, method='lm')['x'][0]
    return pH

def _titrationDefaults():
    """Set titration simulation defaults to resemble Dickson CRM, batch 144."""
    return {
        'acidVolStep': 0.15,
        'alk0': 2238.6e-6,
        'buretteCorrection': 1,
        'concAcid': 0.1,
        'emf0': 660,
        'extraVolAcid': 0,
        'maxVolAcid': 4.1,
        'pSal': 33.571,
        'tempK': 298.15,
        'totalCarbonate': 2031.53e-6,
        'totalPhosphate': 0.31e-6,
        'totalSilicate': 2.5e-6,
        'volSample': 100,
        'totalAmmonia': 0,
        'totalH2Sulfide': 0,
        'WhichKs': 10,
        'WhoseKSO4': 1,
        'WhoseKF': 1,
        'WhoseTB': 2,
    }

def _titrationSettings(**kwargs):
    """Merge user-specified inputs with titration defaults."""
    return {k: v if v is not None else _titrationDefaults()[k]
        for k, v in kwargs.items()}

def titration(acidVolStep=None, alk0=None, buretteCorrection=None,
        concAcid=None, emf0=None, extraVolAcid=None, maxVolAcid=None,
        pSal=None, tempK=None, totalCarbonate=None, totalPhosphate=None,
        totalSilicate=None, volSample=None, totalAmmonia=None,
        totalH2Sulfide=None, WhichKs=None, WhoseKSO4=None, WhoseKF=None,
        WhoseTB=None):
    """Simulate a titration dataset."""
    ts = _titrationSettings(**locals())
    volAcid = (arange(0, ts['maxVolAcid'], ts['acidVolStep']) +
        ts['extraVolAcid'])
    massSample = ts['volSample']*density.sw(ts['tempK'], ts['pSal'])*1e-3
    massAcid = ts['buretteCorrection']*volAcid*density.acid(ts['tempK'])*1e-3
    concTotals = concentrations.concTotals(ts['pSal'], ts['totalCarbonate'],
        ts['totalPhosphate'], ts['totalSilicate'],
        totalAmmonia=ts['totalAmmonia'], totalH2Sulfide=ts['totalH2Sulfide'],
        WhichKs=ts['WhichKs'], WhoseTB=ts['WhoseTB'])
    eqConstants = dissociation.eqConstants(ts['tempK'], ts['pSal'], concTotals,
        WhichKs=ts['WhichKs'], WhoseKSO4=ts['WhoseKSO4'], WhoseKF=ts['WhoseKF'])
    pHSim = pH(massAcid, massSample, ts['concAcid'], ts['alk0'], concTotals,
        eqConstants)
    emf = solve.h2emf(10.0**-pHSim, ts['emf0'], ts['tempK'])
    return volAcid, emf, full_like(emf, ts['tempK'])

def Potentiometric(acidVolStep=None, alk0=None, buretteCorrection=None,
        concAcid=None, emf0=None, extraVolAcid=None, maxVolAcid=None,
        pSal=None, tempK=None, totalCarbonate=None, totalPhosphate=None,
        totalSilicate=None, volSample=None, totalAmmonia=None,
        totalH2Sulfide=None, WhichKs=None, WhoseKSO4=None, WhoseKF=None,
        WhoseTB=None):
    """Simulate a titration dataset as a titration.Potentiometric object."""
    volAcid, emf, tempK = titration(**locals())
    ts = _titrationSettings(**locals())
    return titr.Potentiometric(volAcid, emf, tempK, ts['pSal'], ts['volSample'],
        totalCarbonate=ts['totalCarbonate'],
        totalPhosphate=ts['totalPhosphate'], totalSilicate=ts['totalSilicate'],
        totalAmmonia=ts['totalAmmonia'], totalH2Sulfide=ts['totalH2Sulfide'],
        WhichKs=ts['WhichKs'], WhoseKSO4=ts['WhoseKSO4'], WhoseKF=ts['WhoseKF'],
        WhoseTB=ts['WhoseTB'])
