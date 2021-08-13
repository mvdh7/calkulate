# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019--2021  Matthew P. Humphreys  (GNU GPLv3)
"""Interfaces with external packages."""

import pandas as pd
import PyCO2SYS as pyco2
from . import default


calk_to_pyco2__totals = {
    "total_ammonia": "TNH3",
    "total_phosphate": "TPO4",
    "total_silicate": "TSi",
    "total_sulfide": "TH2S",
    "total_borate": "TB",
    "total_fluoride": "TF",
    "total_sulfate": "TSO4",
    "total_alpha": "total_alpha",
    "total_beta": "total_beta",
}
calk_to_pyco2__k_constants = {
    "k_alpha": "k_alpha",
    "k_ammonia": "KNH3",
    "k_beta": "k_beta",
    "k_bisulfate": "KSO4",
    "k_borate": "KB",
    "k_carbonic_1": "K1",
    "k_carbonic_2": "K2",
    "k_fluoride": "KF",
    "k_phosphoric_1": "KP1",
    "k_phosphoric_2": "KP2",
    "k_phosphoric_3": "KP3",
    "k_silicate": "KSi",
    "k_sulfide": "KH2S",
    "k_water": "KW",
}
pyco2_to_calk__totals = {v: k for k, v in calk_to_pyco2__totals.items()}
pyco2_to_calk__k_constants = {v: k for k, v in calk_to_pyco2__k_constants.items()}


def get_totals(
    salinity,
    dic=0,
    total_alpha=0,
    total_beta=0,
    total_ammonia=0,
    total_phosphate=0,
    total_silicate=0,
    total_sulfide=0,
    total_borate=None,
    total_fluoride=None,
    total_sulfate=None,
    opt_k_carbonic=default.opt_k_carbonic,
    opt_total_borate=default.opt_total_borate,
):
    """Get dict of total substance contents (undiluted) from inputs and PyCO2SYS.

    Inputs in micromol/kg-solution, outputs in mol/kg-solution.
    """
    # Create raw totals dict using PyCO2SYS
    totals_pyco2 = {}
    if total_alpha > 0:
        totals_pyco2["alpha"] = total_alpha * 1e-6
    if total_beta > 0:
        totals_pyco2["beta"] = total_beta * 1e-6
    if not pd.isnull(total_borate):
        totals_pyco2["TB"] = total_borate * 1e-6
    if not pd.isnull(total_fluoride):
        totals_pyco2["TF"] = total_fluoride * 1e-6
    if not pd.isnull(total_sulfate):
        totals_pyco2["TSO4"] = total_sulfate * 1e-6
    totals_pyco2 = pyco2.salts.assemble(
        salinity,
        total_silicate,
        total_phosphate,
        total_ammonia,
        total_sulfide,
        opt_k_carbonic,
        opt_total_borate,
        totals=totals_pyco2,
    )
    # Reorganise totals dict for Calkulate, all in mol/kg-solution
    totals = {}
    totals["dic"] = dic * 1e-6
    if total_alpha > 0:
        totals["total_alpha"] = total_alpha * 1e-6
    if total_beta > 0:
        totals["total_beta"] = total_beta * 1e-6
    if total_ammonia > 0:
        totals["total_ammonia"] = totals_pyco2["TNH3"]
    if total_silicate > 0:
        totals["total_silicate"] = totals_pyco2["TSi"]
    if total_phosphate > 0:
        totals["total_phosphate"] = totals_pyco2["TPO4"]
    if total_sulfide > 0:
        totals["total_sulfide"] = totals_pyco2["TH2S"]
    totals["total_borate"] = totals_pyco2["TB"]
    totals["total_fluoride"] = totals_pyco2["TF"]
    totals["total_sulfate"] = totals_pyco2["TSO4"]
    return totals, totals_pyco2


def get_k_constants(
    totals_pyco2,
    temperature,
    k_alpha=None,
    k_ammonia=None,
    k_beta=None,
    k_bisulfate=None,
    k_borate=None,
    k_carbonic_1=None,
    k_carbonic_2=None,
    k_fluoride=None,
    k_phosphoric_1=None,
    k_phosphoric_2=None,
    k_phosphoric_3=None,
    k_silicate=None,
    k_sulfide=None,
    k_water=None,
    opt_k_bisulfate=default.opt_k_bisulfate,
    opt_k_carbonic=default.opt_k_carbonic,
    opt_k_fluoride=default.opt_k_fluoride,
    opt_pH_scale=default.opt_pH_scale,
    opt_total_borate=default.opt_total_borate,
):
    """Get dict of equilibrium constants from inputs and PyCO2SYS."""
    # Create raw k_constants dict using PyCO2SYS
    k_constants_pyco2 = {"RGas": default.opt_gas_constant * 10}
    if not pd.isnull(k_alpha):
        k_constants_pyco2["k_alpha"] = k_alpha
    if not pd.isnull(k_ammonia):
        k_constants_pyco2["KNH3"] = k_ammonia
    if not pd.isnull(k_beta):
        k_constants_pyco2["k_beta"] = k_beta
    if not pd.isnull(k_bisulfate):
        k_constants_pyco2["KSO4"] = k_bisulfate
    if not pd.isnull(k_borate):
        k_constants_pyco2["KB"] = k_borate
    if not pd.isnull(k_carbonic_1):
        k_constants_pyco2["K1"] = k_carbonic_1
    if not pd.isnull(k_carbonic_2):
        k_constants_pyco2["K2"] = k_carbonic_2
    if not pd.isnull(k_fluoride):
        k_constants_pyco2["KF"] = k_fluoride
    if not pd.isnull(k_phosphoric_1):
        k_constants_pyco2["KP1"] = k_phosphoric_1
    if not pd.isnull(k_phosphoric_2):
        k_constants_pyco2["KP2"] = k_phosphoric_2
    if not pd.isnull(k_phosphoric_3):
        k_constants_pyco2["KP3"] = k_phosphoric_3
    if not pd.isnull(k_silicate):
        k_constants_pyco2["KSi"] = k_silicate
    if not pd.isnull(k_sulfide):
        k_constants_pyco2["KH2S"] = k_sulfide
    if not pd.isnull(k_water):
        k_constants_pyco2["KW"] = k_water
    k_constants_pyco2 = pyco2.equilibria.assemble(
        temperature,
        0,  # assume 1 atm pressure i.e. zero in-water pressure
        totals_pyco2,  # this contains salinity
        opt_pH_scale,
        opt_k_carbonic,
        opt_k_bisulfate,
        opt_k_fluoride,
        None,
        Ks=k_constants_pyco2,
    )
    # Reorganise k_constants dict for Calkulate
    k_constants = {
        k: k_constants_pyco2[v]
        for k, v in calk_to_pyco2__k_constants.items()
        if v in k_constants_pyco2
    }
    # Convert k_fluoride to the Total scale if needed
    if (opt_pH_scale == 1) and "TSO4" in totals_pyco2:
        pH_free_to_total = 1 + totals_pyco2["TSO4"] / k_constants["k_bisulfate"]
        k_constants["k_fluoride"] = k_constants["k_fluoride"] * pH_free_to_total
    return k_constants
