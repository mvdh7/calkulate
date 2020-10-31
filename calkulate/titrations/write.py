import pandas as pd


def to_dat(
    titration,
    filepath,
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
        for _, row in titration.iterrows():
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
                    row.titrant_amount, row.measurement, row.temperature,
                )
            )
