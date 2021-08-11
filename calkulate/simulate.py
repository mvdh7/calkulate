# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019--2021  Matthew P. Humphreys  (GNU GPLv3)
"""Simulate solution properties during a titration."""

import numpy as np
import PyCO2SYS as pyco2
from . import default


def alkalinity_components(pH, totals, k_constants, opt_pH_scale=default.opt_pH_scale):
    """Calculate chemical speciation from pH.

    totals should include dilution correction and be in mol/kg-solution.

    k_constants should be on the scale specified by opt_pH_scale, with:
        1 = Total                pH = -log10([H+] + [HSO4-])
        2 = Seawater             pH = -log10([H+] + [HSO4-] + [HF])
        3 = Free      [default]  pH = -log10([H+])
    Note that when using the Total scale, k_fluoride must also be on the Total scale
    (this does not happen by default in PyCO2SYS!).

    Outputs are substance contents in mol/kg-solution.
    """
    # Check opt_pH_scale is valid
    opt_pH_scales = [1, 2, 3]
    assert (
        opt_pH_scale in opt_pH_scales
    ), "opt_pH_scale must be 1 (Total), 2 (Seawater) or 3 (Free)."
    # Build up dict of solution components
    components = {}
    h = components["H"] = 10.0 ** -pH
    if "k_water" in k_constants:
        components["OH"] = k_constants["k_water"] / h
    if "dic" in totals:
        TCO2 = totals["dic"]
        K1 = k_constants["k_carbonic_1"]
        K2 = k_constants["k_carbonic_2"]
        components["CO2"] = TCO2 / (1 + K1 / h + K1 * K2 / h ** 2)
        components["HCO3"] = K1 * components["CO2"] / h
        components["CO3"] = K2 * components["HCO3"] / h
    if "total_borate" in totals:
        TB = totals["total_borate"]
        KB = k_constants["k_borate"]
        components["BOH4"] = TB * KB / (KB + h)
    if "total_phosphate" in totals:
        TPO4 = totals["total_phosphate"]
        KP1 = k_constants["k_phosphoric_1"]
        KP2 = k_constants["k_phosphoric_2"]
        KP3 = k_constants["k_phosphoric_3"]
        phosphoric_denom = h ** 3 + KP1 * h ** 2 + KP1 * KP2 * h + KP1 * KP2 * KP3
        components["H3PO4"] = TPO4 * h ** 3 / phosphoric_denom
        components["HPO4"] = TPO4 * KP1 * KP2 * h / phosphoric_denom
        components["PO4"] = TPO4 * KP1 * KP2 * KP3 / phosphoric_denom
    if "total_silicate" in totals:
        TSi = totals["total_silicate"]
        KSi = k_constants["k_silicate"]
        components["H3SiO4"] = TSi * KSi / (KSi + h)
    if "total_ammonia" in totals:
        TNH3 = totals["total_ammonia"]
        KNH3 = k_constants["k_ammonia"]
        components["NH3"] = TNH3 * KNH3 / (KNH3 + h)
    if "total_sulfide" in totals:
        TH2S = totals["total_sulfide"]
        KH2S = k_constants["k_sulfide"]
        components["HS"] = TH2S * KH2S / (KH2S + h)
    if "total_alpha" in totals:
        total_alpha = totals["total_alpha"]
        k_alpha = k_constants["k_alpha"]
        alpha = total_alpha * k_alpha / (k_alpha + h)
        alphaH = total_alpha - alpha
        components["alk_alpha"] = np.where(
            -np.log10(k_alpha) <= default.zlp, -alphaH, alpha
        )
    if "total_beta" in totals:
        total_beta = totals["total_beta"]
        k_beta = k_constants["k_beta"]
        beta = total_beta * k_beta / (k_beta + h)
        betaH = total_beta - beta
        components["alk_beta"] = np.where(
            -np.log10(k_beta) <= default.zlp, -betaH, beta
        )
    # pH-scale-dependent components
    if opt_pH_scale in [1, 3]:
        if "total_fluoride" in totals:
            TF = totals["total_fluoride"]
            KF = k_constants["k_fluoride"]
            components["HF"] = TF * h / (KF + h)
    if opt_pH_scale == 3:
        if "total_sulfate" in totals:
            TSO4 = totals["total_sulfate"]
            KSO4 = k_constants["k_bisulfate"]
            components["HSO4"] = TSO4 * h / (KSO4 + h)
    return components


# Multipliers for each component in the alkalinity equation.
# Keys correspond to the output dict of alkalinity_components().
component_multipliers = {
    "H": -1,
    "OH": +1,
    "CO2": 0,
    "HCO3": +1,
    "CO3": +2,
    "BOH4": +1,
    "HSO4": -1,
    "HF": -1,
    "H3PO4": -1,
    "HPO4": +1,
    "PO4": +2,
    "H3SiO4": +1,
    "NH3": +1,
    "HS": +1,
    "alk_alpha": +1,
    "alk_beta": +1,
    "alkalinity_estimate": +1,  # for plotting
}


def alkalinity(pH, totals, k_constants, opt_pH_scale=default.opt_pH_scale):
    """Estimate total alkalinity from pH and total salts in mol/kg-solution."""
    components = alkalinity_components(
        pH, totals, k_constants, opt_pH_scale=opt_pH_scale
    )
    return np.sum([v * component_multipliers[k] for k, v in components.items()], axis=0)
