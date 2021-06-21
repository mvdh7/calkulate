# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019--2021  Matthew P. Humphreys  (GNU GPLv3)
"""Convert between various properties."""

import numpy as np
from . import constants, default, interface

# For convenience
F = constants.faraday
R = constants.ideal_gas


def get_dilution_factor(titrant_mass, analyte_mass):
    """Calculate the factor for dilution of the analyte by the titrant."""
    return analyte_mass / (titrant_mass + analyte_mass)


def dilute_totals(totals, titrant_mass, analyte_mass):
    """Apply the dilution factor to all elements of totals."""
    dilution_factor = get_dilution_factor(titrant_mass, analyte_mass)
    return {k: v * dilution_factor for k, v in totals.items()}


def dilute_totals_pyco2(totals_pyco2, titrant_mass, analyte_mass):
    """Apply the dilution factor to the appropriate elements of totals_pyco2."""
    dilution_factor = get_dilution_factor(titrant_mass, analyte_mass)
    return {
        k: v * dilution_factor
        if k
        in (
            "TB",
            "TF",
            "TSO4",
            "TCa",
            "alpha",
            "beta",
            "TPO4",
            "TSi",
            "TNH3",
            "TH2S",
        )
        else v
        for k, v in totals_pyco2.items()
    }


def totals_to_pyco2(totals, salinity):
    """Convert Calkulate-style totals into PyCO2SYS-style totals_pyco2."""
    totals_pyco2 = {
        p: totals[c] if c in totals else 0.0
        for p, c in interface.pyco2_to_calk__totals.items()
    }
    totals_pyco2["Sal"] = salinity
    return totals_pyco2


def dilute_totals_H2SO4(totals, titrant_molinity, titrant_mass, analyte_mass):
    """Apply the dilution factor to all elements of totals when the titrant is H2SO4."""
    dilution_factor = get_dilution_factor(titrant_mass, analyte_mass)
    totals = {
        k: v * dilution_factor
        if k != "total_sulfate"
        else (v * analyte_mass + titrant_molinity * titrant_mass)
        / (analyte_mass + titrant_mass)
        for k, v in totals.items()
    }
    return totals


def dilute_totals_pyco2_H2SO4(
    totals_pyco2, titrant_molinity, titrant_mass, analyte_mass
):
    """Apply the dilution factor to the appropriate elements of totals_pyco2
    when the titrant is H2SO4.
    """
    dilution_factor = get_dilution_factor(titrant_mass, analyte_mass)
    totals_pyco2 = {
        k: v * dilution_factor
        if k
        in (
            "TB",
            "TF",
            "TCa",
            "alpha",
            "beta",
            "TPO4",
            "TSi",
            "TNH3",
            "TH2S",
        )
        else v
        for k, v in totals_pyco2.items()
    }
    totals_pyco2["TSO4"] = (
        totals_pyco2["TSO4"] * analyte_mass + titrant_molinity * titrant_mass
    ) / (analyte_mass + titrant_mass)
    return totals_pyco2


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


def pH_to_emf(pH, emf0, temperature):
    """Convert pH to EMF."""
    h = 10.0 ** -pH
    return h_to_emf(h, emf0, temperature)


def f_to_demf0(f, temperature):
    """Convert Dickson's f factor into del-EMF0."""
    temperature_K = temperature + constants.absolute_zero
    return np.log(f) * R * temperature_K / F
