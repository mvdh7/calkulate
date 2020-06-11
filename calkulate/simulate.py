# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew P. Humphreys  (GNU GPLv3)
"""Simulate solution properties during a titration."""

import numpy as np


def alkalinity_components(pH, total_salts, equilibrium_constants, total_carbonate=0):
    """Calculate chemical speciation from pH."""
    # Unpack total salts (these should already include the dilution correction)
    TCO2 = total_carbonate
    TB = total_salts["TB"]
    TSO4 = total_salts["TSO4"]
    TF = total_salts["TF"]
    TPO4 = total_salts["TPO4"]
    TSi = total_salts["TSi"]
    TNH3 = total_salts["TNH3"]
    TH2S = total_salts["TH2S"]
    # Unpack equilibrium constants
    KW = equilibrium_constants["KW"]
    K1 = equilibrium_constants["K1"]
    K2 = equilibrium_constants["K2"]
    KB = equilibrium_constants["KB"]
    KSO4 = equilibrium_constants["KSO4"]
    KF = equilibrium_constants["KF"]
    KP1 = equilibrium_constants["KP1"]
    KP2 = equilibrium_constants["KP2"]
    KP3 = equilibrium_constants["KP3"]
    KSi = equilibrium_constants["KSi"]
    KNH3 = equilibrium_constants["KNH3"]
    KH2S = equilibrium_constants["KH2S"]
    # Calculate components
    h = 10.0 ** -pH
    hydroxide = KW / h
    aqueous_CO2 = TCO2 / (1 + K1 / h + K1 * K2 / h ** 2)
    bicarbonate = K1 * aqueous_CO2 / h
    carbonate = K2 * bicarbonate / h
    borate = TB * KB / (KB + h)
    bisulfate = TSO4 * h / (KSO4 + h)
    hydrogen_fluoride = TF * h / (KF + h)
    phosphoric_0 = TPO4 / (1 + KP1 / h + KP1 * KP2 / h ** 2 + KP1 * KP1 * KP3 / h ** 3)
    phosphoric_2 = TPO4 / (h ** 2 / (KP1 * KP2) + h / KP2 + 1 + KP3 / h)
    phosphoric_3 = TPO4 / (
        h ** 3 / (KP1 * KP2 * KP3) + h ** 2 / (KP2 * KP3) + h / KP3 + 1
    )
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


def alkalinity(h, total_salts, equilibrium_constants, total_carbonate=0):
    """Solve for total alkalinity from [H+]."""
    components = alkalinity_components(
        h, total_salts, equilibrium_constants, total_carbonate=total_carbonate
    )
    return np.sum([v * components[k] for k, v in component_multipliers.items()], axis=0)
