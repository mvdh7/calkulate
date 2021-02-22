import pandas as pd


def read_dat(
    filepath_or_buffer,
    titrant_amount_col=0,
    measurement_col=1,
    temperature_col=2,
    skiprows=2,
    header=None,
    **read_dat_kwargs
):
    """Read titration data from a text file."""
    titration = pd.read_table(
        filepath_or_buffer, skiprows=skiprows, header=header, **read_dat_kwargs
    ).rename(
        columns={
            titrant_amount_col: "titrant_amount",
            measurement_col: "measurement",
            temperature_col: "temperature",
        }
    )
    return titration
