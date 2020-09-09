# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew P. Humphreys  (GNU GPLv3)
"""Import and export titration data."""

import numpy as np


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


class DatDict(dict):
    def write(self, fname, **kwargs):
        write_dat(self, fname, **kwargs)


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
