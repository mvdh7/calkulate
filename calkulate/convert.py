# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew P. Humphreys  (GNU GPLv3)
"""Convert between various properties."""

import numpy as np
from . import constants, solve

F = constants.faraday
R = constants.ideal_gas


def dilution_factor(analyte_mass, mixture_mass):
    """Calculate the factor for dilution of the analyte by the titrant."""
    return analyte_mass / mixture_mass


def emf_to_h(emf, emf0, temperature):
    """Convert EMF to [H+]."""
    # DAA03 Eq. (13) with typo corrected (i.e. EMF and EMF0 switched)
    temperature_K = temperature + constants.absolute_zero
    return np.exp((emf - emf0) * F / (R * temperature_K))


def emf_to_pH(emf, emf0, temperature):
    """Convert EMF to pH."""
    return -np.log10(emf_to_h(emf, emf0, temperature))


def h_to_emf(h, emf0, temperature):
    """Convert [H+] to EMF."""
    temperature_K = temperature + constants.absolute_zero
    return emf0 + np.log(h) * R * temperature_K / F


def f_to_dEmf0(f, temperature):
    temperature_K = temperature + constants.absolute_zero
    return np.log(f) * R * temperature_K / F


def measurement_type_to_solver(measurement_type):
    m2s = {
        "EMF": solve.complete_emf,
        "pH": solve.complete_pH,
    }
    return m2s[measurement_type]


calk_to_pyco2 = {
    "total_ammonia": "TNH3",
    "total_phosphate": "TPO4",
    "total_silicate": "TSi",
    "total_sulfide": "TH2S",
    "total_borate": "TB",
    "total_fluoride": "TF",
    "total_sulfate": "TSO4",
    "k_ammonia": "KNH3",
    "k_borate": "KB",
    "k_bisulfate": "KSO4",
    "k_carbonic_1": "K1",
    "k_carbonic_2": "K2",
    "k_fluoride": "KF",
    "k_phosphate_1": "KP1",
    "k_phosphate_2": "KP2",
    "k_phosphate_3": "KP3",
    "k_silicate": "KSi",
    "k_sulfide": "KH2S",
    "k_water": "KW",
}

pyco2_to_calk = {v: k for k, v in calk_to_pyco2.items()}
