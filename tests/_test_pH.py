# %%
import os

import numpy as np
import pandas as pd

import calkulate as calk
from calkulate.dataset import calibrate, calkulate, solve, solve_row


# Create metadata df based on all *.old files from a given path
file_path = "tests/data/ts-tiamo/"
file_name = [f for f in os.listdir(file_path) if f.endswith(".old")]
ds = pd.DataFrame({"file_name": file_name, "file_path": file_path})

# Add extra info about samples
ds["analyte_volume"] = 25  # ml
ds["salinity"] = 33.231
ds["alkalinity_certified"] = 2220.62  # should be NaN if not a CRM

# Calibrate
ds["file_good"] = True
ds["analyte_mass"] = np.nan
ds["dic"] = 2000
ds["measurement_type"] = "pH"
ds["pH_min"] = 2
ds["pH_max"] = 5
ds_row = ds.iloc[0]
prepare_kwargs = {
    "analyte_volume": 25,
    "read_dat_method": "tiamo_de",
}

ds = calk.Dataset(ds)
ds.calkulate(file_type="tiamo_de")

# row = calibrate_row(ds.iloc[0], file_type="tiamo_de")
# calibrate(ds, file_type="tiamo_de")
# row = solve_row(ds.iloc[0], file_type="tiamo_de")
# solve(ds, file_type="tiamo_de")
# calkulate(ds, file_type="tiamo_de")
