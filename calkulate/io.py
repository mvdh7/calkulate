# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew P. Humphreys  (GNU GPLv3)
"""Import and export titration data."""

import numpy as np, pandas as pd
from matplotlib import dates as mdates


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


def add_func_cols(df, func, *args, **kwargs):
    """Add results of df.apply(func) to df as new columns."""
    return df.join(df.apply(lambda x: func(x, *args, **kwargs), axis=1))


def dbs_datetime(dbsx):
    """Convert date and time from .dbs file into datetime and datenum."""
    dspl = dbsx["date"].split("/")
    analysis_datetime = np.datetime64(
        "-".join(("20" + dspl[2], dspl[0], dspl[1])) + "T" + dbsx["time"]
    )
    return pd.Series(
        {
            "analysis_datetime": analysis_datetime,
            "analysis_datenum": mdates.date2num(analysis_datetime),
        }
    )


def get_VINDTA_filenames(dbs):
    """Determine VINDTA filenames assuming defaults were used based on the dbs."""
    dbs["file_name"] = dbs.apply(
        lambda x: "{}-{}  {}  ({}){}.dat".format(
            x.station, x.cast, x.niskin, x.depth, x.bottle
        ),
        axis=1,
    )
    return dbs


def read_dbs(fname, analyte_volume=100.0, analyte_mass=None, file_path=None):
    """Import one .dbs file as single DataFrame."""
    headers = np.genfromtxt(fname, delimiter="\t", dtype=str, max_rows=1)
    dbs = pd.read_table(fname, header=0, names=headers, usecols=headers)
    dbs["dbs_fname"] = fname
    dbs = add_func_cols(dbs, dbs_datetime)
    dbs = get_VINDTA_filenames(dbs)
    if analyte_mass is None:
        dbs["analyte_volume"] = analyte_volume
    else:
        dbs["analyte_mass"] = analyte_mass
    if file_path is not None:
        dbs["file_path"] = file_path
    return dbs
