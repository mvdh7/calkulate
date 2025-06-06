# %%
import os

import pandas as pd

import calkulate as calk


def test_tiamo_ds():
    # Create the metadata dataframe (could make this another way,
    # e.g. in a spreadsheet that is then imported with pandas or calkulate)
    ds = {}
    ds["file_path"] = "tests/data/ts-tiamo/"
    ds["file_name"] = [
        f for f in os.listdir(ds["file_path"]) if f.endswith(".old")
    ]
    ds["analyte_volume"] = 25  # ml
    ds["salinity"] = 33.231
    ds["total_silicate"] = 5.7
    ds["total_phosphate"] = 0.56
    ds["dic"] = 2046.37
    ds["alkalinity_certified"] = 2220.62
    ds = pd.DataFrame(ds)
    # Calibrate and solve
    # kwargs_tiamo could alternatively be added as columns to the ds
    kwargs_tiamo = dict(
        titrant_molinity_init=0.01,
        file_type="tiamo_de",
        solve_mode="pH_gran",  # can try "pH" or "pH_adjust" or "pH_gran"
    )
    calk.calibrate(ds, **kwargs_tiamo)
    assert ds.alkalinity.notnull().all()


# test_tiamo_ds()
