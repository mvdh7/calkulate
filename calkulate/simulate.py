# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew P. Humphreys  (GNU GPLv3)
"""Simulate solution properties during a titration."""

import numpy as np


def alkalinity_components(pH, totals, k_constants, dic=0):
    """Calculate chemical speciation from pH."""
    # Unpack total salts (these should already include the dilution correction)
    TCO2 = dic
    TB = totals["TB"]
    TSO4 = totals["TSO4"]
    TF = totals["TF"]
    TPO4 = totals["TPO4"]
    TSi = totals["TSi"]
    TNH3 = totals["TNH3"]
    TH2S = totals["TH2S"]
    # Unpack equilibrium constants (should be on the Free scale)
    KW = k_constants["KW"]
    K1 = k_constants["K1"]
    K2 = k_constants["K2"]
    KB = k_constants["KB"]
    KSO4 = k_constants["KSO4"]
    KF = k_constants["KF"]
    KP1 = k_constants["KP1"]
    KP2 = k_constants["KP2"]
    KP3 = k_constants["KP3"]
    KSi = k_constants["KSi"]
    KNH3 = k_constants["KNH3"]
    KH2S = k_constants["KH2S"]
    # Calculate components
    h = 10.0 ** -pH
    hydroxide = KW / h
    aqueous_CO2 = TCO2 / (1 + K1 / h + K1 * K2 / h ** 2)
    bicarbonate = K1 * aqueous_CO2 / h
    carbonate = K2 * bicarbonate / h
    borate = TB * KB / (KB + h)
    bisulfate = TSO4 * h / (KSO4 + h)
    hydrogen_fluoride = TF * h / (KF + h)
    phosphoric_denom = h ** 3 + KP1 * h ** 2 + KP1 * KP2 * h + KP1 * KP2 * KP3
    phosphoric_0 = TPO4 * h ** 3 / phosphoric_denom
    phosphoric_2 = TPO4 * KP1 * KP2 * h / phosphoric_denom
    phosphoric_3 = TPO4 * KP1 * KP2 * KP3 / phosphoric_denom
    silicate = TSi * KSi / (KSi + h)
    ammonia = TNH3 * KNH3 / (KNH3 + h)
    hydrogen_sulfide = TH2S * KH2S / (KH2S + h)
    return {
        "H": h,
        "OH": hydroxide,
        "CO2aq": aqueous_CO2,
        "HCO3": bicarbonate,
        "CO3": carbonate,
        "BOH4": borate,
        "HSO4": bisulfate,
        "HF": hydrogen_fluoride,
        "H3PO4": phosphoric_0,
        "HPO4": phosphoric_2,
        "PO4": phosphoric_3,
        "SiOOH3": silicate,
        "NH3": ammonia,
        "HS": hydrogen_sulfide,
    }


# Multipliers for each component in the alkalinity equation.
# Keys correspond to the output dict of alkalinity_components().
component_multipliers = {
    "H": -1,
    "OH": +1,
    "HCO3": +1,
    "CO3": +2,
    "BOH4": +1,
    "HSO4": -1,
    "HF": -1,
    "H3PO4": -1,
    "HPO4": +1,
    "PO4": +2,
    "SiOOH3": +1,
    "NH3": +1,
    "HS": +1,
}


def alkalinity(pH, totals, k_constants, dic=0):
    """Estimate total alkalinity from [H+] and total salts."""
    components = alkalinity_components(pH, totals, k_constants, dic=dic)
    return np.sum([v * components[k] for k, v in component_multipliers.items()], axis=0)
