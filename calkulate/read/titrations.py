# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019--2025  Matthew P. Humphreys  (GNU GPLv3)
"""Import data for Calkulate and export its results."""

import re
from collections import namedtuple
from warnings import warn

import numpy as np
import pandas as pd


DatData = namedtuple(
    "DatData", ("titrant_amount", "measurement", "temperature")
)


def read_dat_genfromtxt(
    file_name,
    col_titrant_amount=0,
    col_measurement=1,
    col_temperature=2,
    delimiter="\t",
    skip_header=2,
    **kwargs,
):
    """Import a titration dataset from a .dat file with numpy.genfromtxt."""
    data = np.genfromtxt(
        file_name, delimiter=delimiter, skip_header=skip_header, **kwargs
    )
    titrant_amount = data[:, col_titrant_amount]
    measurement = data[:, col_measurement]
    temperature = data[:, col_temperature]
    return DatData(titrant_amount, measurement, temperature)


def read_dat_pclims(
    file_name,
    col_titrant_amount=1,
    col_measurement=2,
    col_temperature=5,
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
    titrant_amount = data[:, col_titrant_amount]
    measurement = data[:, col_measurement]
    temperature = data[:, col_temperature]
    return DatData(titrant_amount, measurement, temperature)


def read_dat_orgalk_excel(file_name):
    """Import a titration dataset from an Excel file formatted for the NIOZ
    organic alkalinity project.

    Parameters
    ----------
    file_name : str
        The file name (and path).

    Returns
    -------
    acid_amount : array-like
        The amount of acid added at each titration step (should be in ml).
    base_amount : array-like
        The amount of base added at each titration step (should be in ml).
    emf : array-like
        The EMF at each titration step (should be in mV).
    temperature : array-like
        The temperature at each titration step (should be in °C).
    """
    data = pd.read_excel(file_name, skiprows=1)
    return (
        data.acid.to_numpy(),
        data.base.to_numpy(),
        data.emf.to_numpy(),
        data.temperature.to_numpy(),
    )


def _read_tiamo_de_df(file_name, encoding="unicode_escape"):
    with open(file_name, "r", encoding=encoding) as f:
        data = f.read().splitlines()
    gran_line = data.index("Gran.1")
    volume_start = float(data[gran_line - 3].split(";")[2])
    data = pd.read_csv(
        file_name,
        skiprows=gran_line + 1,
        encoding=encoding,
        sep=";",
    ).rename(
        columns={
            "Volumen [mL]": "volume",
            "Messwert": "pH",
            "dMW [pH]": "dpH",
            "Zeit [s]": "time",
            "Temperatur [°C]": "temperature",
            "Berechnet 1": "b1",
            "Berechnet 2": "b2",
            "Berechnet 3": "b3",
        }
    )
    data["volume"] += volume_start
    return data


def read_tiamo_de(file_name, encoding="unicode_escape"):
    """Read volume, pH and temperature from a German Tiamo data file
    (as received from T. Steinhoff, GEOMAR).

    Parameters
    ----------
    file_name : str
        The file name (and path).
    encoding : str, optional
        The file encoding, by default "unicode_escape"

    Returns
    -------
    TiamoData
        A named tuple containing the fields
            volume : float
                Volume of titrant added in ml.
            pH : float
                pH measured by the electrode.
            temperature : float
                Temperature in °C.
    """
    data = _read_tiamo_de_df(file_name, encoding=encoding)
    return DatData(data.volume.values, data.pH.values, data.temperature.values)


file_types = {
    "genfromtxt": read_dat_genfromtxt,
    "vindta": read_dat_genfromtxt,
    "orgalk_excel": read_dat_orgalk_excel,
    "pclims": read_dat_pclims,
    "tiamo_de": read_tiamo_de,
}
keys_read_dat = [
    "col_measurement",
    "col_temperature",
    "col_titrant_amount",
    "delimiter",
    "encoding",
    "file_type",
    "n_cols",
    "skip_header",
]


def read_dat(file_name, file_type="genfromtxt", **kwargs):
    """Import a titration dataset from a .dat file."""
    if file_type not in file_types:
        file_type = "genfromtxt"
        warn(f"method '{file_type}' not recognised, trying 'genfromtxt'.")
    dat_data = file_types[file_type](file_name, **kwargs)
    return dat_data


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
    **open_kwargs,
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
