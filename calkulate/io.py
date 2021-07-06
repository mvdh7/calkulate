# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019--2021  Matthew P. Humphreys  (GNU GPLv3)
"""Import data for Calkulate and export its results."""

import re
import numpy as np, pandas as pd
from matplotlib import dates as mdates
from . import dataset, default


def read_dat_genfromtxt(
    file_name,
    titrant_amount_col=0,
    measurement_col=1,
    temperature_col=2,
    delimiter="\t",
    skip_header=2,
    **kwargs
):
    """Import a titration dataset from a .dat file with numpy.genfromtxt."""
    data = np.genfromtxt(
        file_name, delimiter=delimiter, skip_header=skip_header, **kwargs
    )
    titrant_amount = data[:, titrant_amount_col]
    measurement = data[:, measurement_col]
    temperature = data[:, temperature_col]
    return titrant_amount, measurement, temperature


def read_dat_pclims(
    file_name,
    titrant_amount_col=1,
    measurement_col=2,
    temperature_col=5,
    n_cols=6,
):
    """Import a titration dataset from a PC LIMS Report."""
    # Import full data file
    with open(file_name, mode="r", encoding="unicode-escape") as f:
        data_full = f.read().splitlines()
    # Select only the lines containing the titration data
    re_data = re.compile(r"[\d\.]*" + r"\t[\d\.]*" * (n_cols - 1))
    data = []
    i = 0
    stop = False
    found_any = False
    while not stop and i < len(data_full):
        line = data_full[i]
        if re_data.match(line):
            found_any = True
            data.append(line.split("\t"))
        else:
            if found_any:
                stop = True
        i += 1
    # Convert data to float arrays
    data = np.array(data, dtype="float64")
    titrant_amount = data[:, titrant_amount_col]
    measurement = data[:, measurement_col]
    temperature = data[:, temperature_col]
    return titrant_amount, measurement, temperature


methods = {"genfromtxt": read_dat_genfromtxt, "pclims": read_dat_pclims}


def read_dat(file_name, method=default.read_dat_method, **kwargs):
    """Import a titration dataset from a .dat file."""
    if method not in methods:
        method = "genfromtxt"
        print("method '{}' not recognised, using 'genfromtxt'.".format(method))
    titrant_amount, measurement, temperature = methods[method](file_name, **kwargs)
    return titrant_amount, measurement, temperature


def write_dat(
    filepath,
    titrant_amount,
    measurement,
    temperature,
    line0="Titration data exported by Calkulate",
    line1="titrant_amount\tmeasurement\ttemperature",
    mode="x",
    titrant_amount_fmt=".3f",
    measurement_fmt=".3f",
    temperature_fmt=".3f",
    **open_kwargs
):
    """Write titration data to a text file."""
    with open(filepath, mode=mode, **open_kwargs) as f:
        f.write("{}\n{}\n".format(line0, line1))
        for titrant_amount_i, measurement_i, temperature_i in zip(
            titrant_amount, measurement, temperature
        ):
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
                    titrant_amount_i,
                    measurement_i,
                    temperature_i,
                )
            )


def read_clipboard(**read_clipboard_kwargs):
    """Import a metadata table from clipboard and pass to read_csv."""
    return dataset.Dataset(pd.read_clipboard(**read_clipboard_kwargs))


def read_csv(filepath_or_buffer, **read_csv_kwargs):
    """Import a metadata table from a CSV file."""
    return dataset.Dataset(pd.read_csv(filepath_or_buffer, **read_csv_kwargs))


def read_excel(io, **read_excel_kwargs):
    """Import a metadata table from an Excel file."""
    return dataset.Dataset(pd.read_excel(io, **read_excel_kwargs))


def read_fwf(filepath_or_buffer, **read_fwf_kwargs):
    """Import a metadata table from fixed-width formatted lines."""
    return dataset.Dataset(pd.read_fwf(filepath_or_buffer, **read_fwf_kwargs))


def read_table(filepath_or_buffer, **read_csv_kwargs):
    """Import a metadata table from a general delimited file."""
    return dataset.Dataset(pd.read_table(filepath_or_buffer, **read_csv_kwargs))


def add_func_cols(df, func, *args, **kwargs):
    """Add results of df.apply(func) to df as new columns."""
    return df.join(df.apply(lambda x: func(x, *args, **kwargs), axis=1))


def dbs_datetime(dbs_row):
    """Convert date and time from .dbs file into datetime."""
    try:
        dspl = dbs_row["date"].split("/")
        analysis_datetime = np.datetime64(
            "-".join(("20" + dspl[2], dspl[0], dspl[1])) + "T" + dbs_row["time"]
        )
    except AttributeError:
        analysis_datetime = np.datetime64("NaT")
    return pd.Series(
        {
            "analysis_datetime": analysis_datetime,
        }
    )


def get_VINDTA_filenames(dbs):
    """Determine VINDTA filenames, assuming defaults were used, based on the dbs."""
    dbs["file_name"] = dbs.apply(
        lambda x: "{}-{}  {}  ({}){}.dat".format(
            int(x.station), int(x.cast), int(x.niskin), x.depth, x.bottle
        ),
        axis=1,
    )
    return dbs


def read_dbs(fname, analyte_volume=100.0, analyte_mass=None, file_path=None):
    """Import one .dbs file from a VINDTA as single DataFrame."""
    headers = np.genfromtxt(fname, delimiter="\t", dtype=str, max_rows=1)
    dbs = pd.read_table(fname, header=0, names=headers, usecols=headers)
    dbs["dbs_fname"] = fname
    dbs = add_func_cols(dbs, dbs_datetime)
    dbs["analysis_datenum"] = mdates.date2num(dbs.analysis_datetime)
    dbs = get_VINDTA_filenames(dbs)
    if analyte_mass is None:
        dbs["analyte_volume"] = analyte_volume
    else:
        dbs["analyte_mass"] = analyte_mass
    if file_path is not None:
        assert isinstance(file_path, str), "file_path must be a string."
        dbs["file_path"] = file_path
    return dataset.Dataset(dbs)
