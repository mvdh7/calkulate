import pandas as pd, numpy as np
from . import Titration


def read_dat(
    filepath_or_buffer,
    titrant_amount_col=0,
    measurement_col=1,
    temperature_col=2,
    skiprows=2,
    header=None,
    method="pandas",
    fread_start_line="$S Mode 1\t01\tDET U\tV1.0",
    **read_dat_kwargs
):
    """Read titration data from a text file."""
    assert method in ["pandas", "fread"], "Method must be 'pandas' or 'fread'."
    if method == "pandas":
        titration = Titration(
            pd.read_table(
                filepath_or_buffer, skiprows=skiprows, header=header, **read_dat_kwargs
            ).rename(
                columns={
                    titrant_amount_col: "titrant_amount",
                    measurement_col: "measurement",
                    temperature_col: "temperature",
                }
            )
        )
    elif method == "fread":
        with open(filepath_or_buffer, mode="r", encoding="ISO-8859-1") as f:
            raw_text = f.read().splitlines()
        ix = np.where(np.array(raw_text) == fread_start_line)[0][0] + 1
        data = []
        line = raw_text[ix]
        while line != "$E":
            data.append(line.split("\t"))
            ix += 1
            line = raw_text[ix]
        data = np.array(data, dtype=float)
        titration = Titration(
            data={
                "titrant_amount": data[:, titrant_amount_col],
                "measurement": data[:, measurement_col],
                "temperature": data[:, temperature_col],
            }
        )
    return titration


kwargs_TiTouch = dict(
    method="fread", titrant_amount_col=1, measurement_col=2, temperature_col=5,
)
