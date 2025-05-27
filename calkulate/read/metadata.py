import matplotlib.dates as mdates
import numpy as np
import pandas as pd

from ..classes import Dataset


def read_clipboard(**read_clipboard_kwargs):
    """Import a metadata table from clipboard and pass to read_csv."""
    return Dataset(pd.read_clipboard(**read_clipboard_kwargs))


def read_csv(filepath_or_buffer, **read_csv_kwargs):
    """Import a metadata table from a CSV file."""
    return Dataset(pd.read_csv(filepath_or_buffer, **read_csv_kwargs))


def read_excel(io, **read_excel_kwargs):
    """Import a metadata table from an Excel file."""
    return Dataset(pd.read_excel(io, **read_excel_kwargs))


def read_fwf(filepath_or_buffer, **read_fwf_kwargs):
    """Import a metadata table from fixed-width formatted lines."""
    return Dataset(pd.read_fwf(filepath_or_buffer, **read_fwf_kwargs))


def read_table(filepath_or_buffer, **read_csv_kwargs):
    """Import a metadata table from a general delimited file."""
    return Dataset(pd.read_table(filepath_or_buffer, **read_csv_kwargs))


def add_func_cols(df, func, *args, **kwargs):
    """Add results of df.apply(func) to df as new columns."""
    return df.join(df.apply(lambda x: func(x, *args, **kwargs), axis=1))


def dbs_datetime(dbs_row):
    """Convert date and time from .dbs file into datetime."""
    try:
        dspl = dbs_row["date"].split("/")
        analysis_datetime = np.datetime64(
            "-".join(("20" + dspl[2], dspl[0], dspl[1]))
            + "T"
            + dbs_row["time"]
        )
    except AttributeError:
        analysis_datetime = np.datetime64("NaT")
    return pd.Series(
        {
            "analysis_datetime": analysis_datetime,
        }
    )


def get_VINDTA_filenames(dbs):
    """Determine VINDTA filenames, assuming defaults were used, based on the
    dbs.
    """
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
    return Dataset(dbs)
