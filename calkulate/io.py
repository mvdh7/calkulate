# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew P. Humphreys  (GNU GPLv3)
"""Import and export titration data."""

import numpy as np


def read_dat(
    fname,
    titrant_amount_col=0,
    measurement_col=1,
    temperature_col=2,
    delimiter="\t",
    skip_header=2,
    **kwargs
):
    """Import titration data from a text file."""
    data = np.genfromtxt(fname, delimiter=delimiter, skip_header=skip_header, **kwargs)
    return {
        "titrant_amount": data[:, titrant_amount_col],
        "mixture_measurement": data[:, measurement_col],
        "mixture_temperature": data[:, temperature_col],
    }


def write_dat(
    titration,
    fname,
    line0="",
    line1="titrant_amount\tmeasurement\ttemperature",
    mode="x",
    titrant_amount_fmt=".5f",
    measurement_fmt=".5f",
    temperature_fmt=".3f",
    **kwargs
):
    """Write titration data to a text file."""
    if titration.measurement_type == "EMF":
        measurement = titration.mixture.emf
    elif titration.measurement_type == "pH":
        measurement = titration.mixture.pH
    if hasattr(titration.titrant, "volume"):
        titrant_amount = titration.titrant.volume
    else:
        titrant_amount = titration.titrant.mass
    with open(fname, mode=mode, **kwargs) as f:
        f.write("{}\n{}\n".format(line0, line1))
        for i in range(len(titration.titrant.volume)):
            f.write(
                (
                    "{:"
                    + titrant_amount_fmt
                    + "}\t{:"
                    + measurement_fmt
                    + "}\t{:"
                    + temperature_fmt
                    + "}\n"
                ).format(
                    titrant_amount[i], measurement[i], titration.mixture.temperature[i],
                )
            )


def check_set(collection, field, default):
    """Check if a field is in a collection and that it is not a NaN.
    Return the field if so, or the default if not.
    """
    if field in collection:
        cf = collection[field]
        if isinstance(cf, str):
            if cf == "":
                return default
            else:
                return cf
        elif ~np.isnan(cf):
            return cf
        else:
            return default
    else:
        return default
