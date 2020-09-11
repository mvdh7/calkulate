import pandas as pd
from . import Titration


def read_dat(
    filepath_or_buffer,
    titrant_amount_col=0,
    measurement_col=1,
    temperature_col=2,
    skiprows=2,
    header=None,
    **read_table_kwargs
):
    """Read titration data from a text file."""
    return Titration(
        pd.read_table(
            filepath_or_buffer, skiprows=skiprows, header=header, **read_table_kwargs
        ).rename(
            columns={
                titrant_amount_col: "titrant_amount",
                measurement_col: "measurement",
                temperature_col: "temperature",
            }
        )
    )
