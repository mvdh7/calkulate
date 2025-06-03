# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019--2025  Matthew P. Humphreys  (GNU GPLv3)
"""Convert between various properties."""

from collections import namedtuple

import numpy as np
import pandas as pd

from . import constants, density, interface
from .meta import _get_kwarg_keys


Converted = namedtuple(
    "ConvertResult",
    ("titrant_mass", "measurement", "temperature", "analyte_mass", "salinity"),
)

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
        k: (
            v * dilution_factor
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
        )
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
        k: (
            v * dilution_factor
            if k != "total_sulfate"
            else (v * analyte_mass + titrant_molinity * titrant_mass)
            / (analyte_mass + titrant_mass)
        )
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
        k: (
            v * dilution_factor
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
        )
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
    h = 10.0**-pH
    return h_to_emf(h, emf0, temperature)


def f_to_demf0(f, temperature):
    """Convert Dickson's f factor into del-EMF0."""
    temperature_K = temperature + constants.absolute_zero
    return np.log(f) * R * temperature_K / F


def amount_units(
    dat_data,
    salinity,
    analyte_mass=None,
    analyte_volume=None,
    molinity_H2SO4=0.1,
    molinity_HCl=0.1,
    molinity_NaCl=0.6,
    temperature_override=None,
    titrant_amount_unit="ml",
    titrant_density=None,
    titrant="HCl",
):
    """Convert titrant and analyte units to mass in kg.

    Parameters
    ----------
    dat_data : DatData
        Data from a titration file imported with `read.read_dat`.
    salinity : float
        Practical salinity of the analyte.
    analyte_mass : float, optional
        Analyte mass in kg.  Either this or `analyte_volume` must be given.
    analyte_volume : float, optional
        Analyte volume in ml, which is converted to kg assuming analyte is
        seawater.  Either this or `analyte_mass` must be given.
    molinity_H2SO4 : float, optional
        H2SO4 titrant molinity in mol/kg-sol for density calculation, by
        default 0.1.
    molinity_HCl : float, optional
        HCl titrant molinity in mol/kg-sol for density calculation, by
        default 0.1.
    molinity_NaCl : float, optional
        NaCl molinity in mol/kg-sol in the titrant where it is an HCl-NaCl
        mixture, for density calculation, by default 0.6.
    temperature_override : float, optional
        A temperature in °C to use instead of the temperature data in the
        titration file, by default `None`.
    titrant_amount_unit : str, optional
        Units for the titrant amount in the file, one of "ml" (default), "g" or
        "kg".
    titrant_density : float, optional
        Titrant density in g/ml, by default `None`, in which case it is
        calculated from the molinities provided.
    titrant : str, optional
        Which titrant was used, "HCl" (default) or "H2SO4".

    Returns
    -------
    Converted - namedtuple containing the fields
        titrant_mass : array-like float
            Titrant mass through the titration in kg.
        measurement : array-like float
            EMF through the titration in mV, or pH.
        temperature : array-like float
            Temperature through the titration in °C.
        analyte_mass : float
            Mass of the analyte in kg.
    """
    dd = dat_data
    # Overwrite temperature, if requested
    if temperature_override is not None:
        temperature = np.full_like(dd.temperature, temperature_override)
    else:
        temperature = dd.temperature
    # Get titrant mass
    assert titrant_amount_unit.lower() in ["ml", "g", "kg"]
    if titrant_amount_unit.lower() == "ml":
        if titrant_density is None:
            assert titrant.upper() in ["H2SO4", "HCL"]
            if titrant.upper() == "H2SO4":
                titrant_mass = (
                    dd.titrant_amount
                    * density.H2SO4_25C_EAIM(molinity_H2SO4)
                    * 1e-3
                )
            else:
                titrant_mass = (
                    dd.titrant_amount
                    * density.HCl_NaCl_25C_DSC07(
                        molinity_HCl=molinity_HCl,
                        molinity_NaCl=molinity_NaCl,
                    )
                    * 1e-3
                )
        else:
            titrant_mass = dd.titrant_amount * titrant_density * 1e-3
    elif titrant_amount_unit.lower() == "g":
        titrant_mass = dd.titrant_amount * 1e-3
    elif titrant_amount_unit.lower() == "kg":
        titrant_mass = dd.titrant_amount
    # Convert analyte_mass to analyte_volume if necessary
    if pd.isnull(analyte_mass):
        assert not pd.isnull(analyte_volume), (
            "Either `analyte_mass` or `analyte_volume` must be provided"
        )
        analyte_mass = (
            analyte_volume
            * density.seawater_1atm_MP81(
                temperature=temperature[0],
                salinity=salinity,
            )
            * 1e-3
        )
    return Converted(
        titrant_mass, dd.measurement, temperature, analyte_mass, salinity
    )


keys_cau = _get_kwarg_keys(amount_units)
