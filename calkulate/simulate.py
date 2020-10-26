# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew P. Humphreys  (GNU GPLv3)
"""Simulate solution properties during a titration."""

import numpy as np


def alkalinity_components(pH, titration):
    """Calculate chemical speciation from pH."""
    # Unpack total salts (these should already include any dilution correction)
    TCO2 = titration["dic"] * 1e-6
    TB = titration["total_borate"] * 1e-6
    TSO4 = titration["total_sulfate"] * 1e-6
    TF = titration["total_fluoride"] * 1e-6
    TPO4 = titration["total_phosphate"] * 1e-6
    TSi = titration["total_silicate"] * 1e-6
    TNH3 = titration["total_ammonia"] * 1e-6
    TH2S = titration["total_sulfide"] * 1e-6
    # Unpack equilibrium constants (should all be on the Free scale)
    KW = titration["k_water"]
    K1 = titration["k_carbonic_1"]
    K2 = titration["k_carbonic_2"]
    KB = titration["k_borate"]
    KSO4 = titration["k_bisulfate"]
    KF = titration["k_fluoride"]
    KP1 = titration["k_phosphoric_1"]
    KP2 = titration["k_phosphoric_2"]
    KP3 = titration["k_phosphoric_3"]
    KSi = titration["k_silicate"]
    KNH3 = titration["k_ammonia"]
    KH2S = titration["k_sulfide"]
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
        "CO2": aqueous_CO2,
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


def alkalinity(pH, titration):
    """Estimate total alkalinity from [H+] and total salts."""
    components = alkalinity_components(pH, titration)
    return np.sum([v * components[k] for k, v in component_multipliers.items()], axis=0)
