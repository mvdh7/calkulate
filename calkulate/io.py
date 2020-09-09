# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew P. Humphreys  (GNU GPLv3)
"""Import and export titration data."""

import numpy as np, pandas as pd
from matplotlib import dates as mdates


class DatDict(dict):
    def write(self, fname, **kwargs):
        write_dat(self, fname, **kwargs)


class Dataset(pd.DataFrame):
    def get_dat_files(self, **read_dat_kwargs):
        """(Re-)import all .dat files."""
        if "file_good" not in self:
            self["file_good"] = True
        if "dat_dict" in self:
            self.drop(columns="dat_dict", inplace=True)
        dats = {}
        for i, row in self.iterrows():
            if row.file_good:
                if "file_path" in row:
                    fname = row.file_path + row.file_name
                else:
                    fname = row.file_name
                try:
                    dats[i] = read_dat(fname, **read_dat_kwargs)
                except IOError:
                    print("Can't find file: '{}'.".format(fname))
                    dats[i] = None
                except:
                    print("Error importing file: '{}'.".format(fname))
                    dats[i] = None
        return Dataset(self.join(pd.DataFrame({"dat_dict": dats})))


def write_dat(
    dat_dict,
    fname,
    line0="Titration data exported by Calkulate",
    line1="titrant_amount\tmixture_measurement\tmixture_temperature",
    mode="x",
    titrant_amount_fmt=".5f",
    measurement_fmt=".5f",
    temperature_fmt=".3f",
    **kwargs
):
    """Write titration data to a text file."""
    with open(fname, mode=mode, **kwargs) as f:
        f.write("{}\n{}\n".format(line0, line1))
        for i in range(dat_dict["titrant_amount"].size):
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
                    dat_dict["titrant_amount"][i],
                    dat_dict["mixture_measurement"][i],
                    dat_dict["mixture_temperature"][i],
                )
            )


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
    dat_dict = DatDict(
        {
            "titrant_amount": data[:, titrant_amount_col],
            "mixture_measurement": data[:, measurement_col],
            "mixture_temperature": data[:, temperature_col],
        }
    )
    return dat_dict


def read_csv(filepath_or_buffer, **kwargs):
    """Import a Dataset from an Excel file using pandas.read_csv."""
    return Dataset(pd.read_csv(filepath_or_buffer, **kwargs))


def read_excel(*args, **kwargs):
    """Import a Dataset from an Excel file using pandas.read_excel."""
    return Dataset(pd.read_excel(*args, **kwargs))


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
    """Determine VINDTA filenames, assuming defaults were used, based on the dbs."""
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
        assert isinstance(file_path, str), "file_path must be a string."
        dbs["file_path"] = file_path
    return Dataset(dbs)
