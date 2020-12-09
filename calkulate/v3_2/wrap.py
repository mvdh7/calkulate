# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew P. Humphreys  (GNU GPLv3)
import numpy as np
from . import density, io, solve


def _get_from_metadata(metadata, key, default):
    if key in metadata:
        return metadata[key]
    else:
        return default


def get_titrant_mass(
    titration,
    titrant_amount_unit,
    density_func=density.HCl_NaCl_25C_DSC07,
    **density_func_kwargs
):
    """Get mass of titrant through a titration in kg."""
    assert titrant_amount_unit in [
        "ml",
        "g",
    ], "titrant_amount_unit must be 'ml' or 'g'."
    if titrant_amount_unit == "g":
        titration["titrant_mass"] = titration["titrant_amount"] * 1e-3
    elif titrant_amount_unit == "ml":
        titration["titrant_mass"] = (
            titration["titrant_amount"] * density_func(**density_func_kwargs) * 1e-3
        )
    return titration


def prepare_metadata(metadata):
    """Prepare metadata for analysis."""
    assert "file_name" in metadata
    assert "salinity" in metadata
    # Get analyte mass in kg
    if "analyte_mass" not in metadata:
        assert "analyte_volume" in metadata
        metadata["analyte_mass"] = np.nan
    if "analyte_volume" in metadata:
        metadata["analyte_mass"] = np.where(
            np.isnan(metadata["analyte_mass"]),
            metadata["analyte_volume"]
            * density.seawater_1atm_MP81(
                salinity=metadata["salinity"], temperature=25
            ),  # fix temperature!
            metadata["analyte_mass"],
        )
    else:
        assert "analyte_mass" in metadata
    return metadata


def prepare_titration(metadata, titration=None):
    """Prepare titration data for analysis."""
    # Import titration data, if required
    if titration is None:
        titration = io.read_dat(metadata["file_name"])
    # Get values from metadata or defaults if not provided
    titrant_amount_unit = _get_from_metadata(metadata, "titrant_amount_unit", "ml")
    measurement_type = _get_from_metadata(metadata, "measurement_type", "emf")
    titrant_density_func = _get_from_metadata(
        metadata, "titrant_density_func", density.HCl_NaCl_25C_DSC07
    )
    titrant_density_func_kwargs = _get_from_metadata(
        metadata, "titrant_density_func_kwargs", {}
    )
    for k in ["molinity_HCl", "molinity_NaCl"]:
        if k in metadata:
            titrant_density_func_kwargs[k] = metadata[k]
    # Calculate titrant and mixture masses
    titration = get_titrant_mass(
        titration,
        titrant_amount_unit,
        density_func=titrant_density_func,
        **titrant_density_func_kwargs
    )
    titration["mixture_mass"] = titration["titrant_mass"] + metadata["analyte_mass"]
    # Assign measurement type
    assert measurement_type.lower() in [
        "emf",
        "ph",
    ], "measurement_type must be 'emf' of 'pH'."
    if measurement_type.lower() == "emf":
        titration["emf"] = titration["measurement"]
    elif measurement_type.lower() == "ph":
        titration["pH"] = titration["measurement"]
    return titration


def get_gran_estimator(metadata, titration=None, **kwargs):
    """Simple Gran-plot estimator following DAA03 eq. 10."""
    titration = prepare_titration(metadata, titration)
    return solve.gran_estimator(
        titration["mixture_mass"], titration["emf"], titration["temperature"], **kwargs
    )
