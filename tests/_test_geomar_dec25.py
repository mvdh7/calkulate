# %%
import os

import pandas as pd

import calkulate as calk


# Set the acid concentration (mol/kg) to use for "uncalibrated" alkalinity
titrant_molinity_guess = 0.01

# Set the file paths where the data files can be found
path_to_files = "tests/data/geomar-dec25/"
metadata_file = "Acidbatch 4.xlsx"
subpath = "Acidbatch 4 20240712 to 20240725"
#
# This assumes a folder structure like the following:
# .
# └── tests/data/geomar-dec25/              (<-- defined by `path_to_files`)
#     ├── Acidbatch 4.xlsx                  (<-- defined by `metadata_file`)
#     └── Acidbatch 4 20240712 to 20240725/ (<-- defined by `subpath`)
#         ├── 20240712/
#         │   ├── Junk-SW-01.old
#         │   ├── Junk-SW-02.old
#         │   └── ...
#         ├── 20240716/
#         │   └── ...
#         └── ...
#
# The subfolder names (20240712, 20240716, etc.) need to be specified in the
# `session` column of the metadata spreadsheet.
#
# The file names, excluding the ".old", need to be in the `name` column.
#
#            ==================================================
#
#             You shouldn't need to change anything below here
#                 (but there is more info at the bottom)
#
#            ==================================================
#
# Next, we use the information provided above to construct columns containing
# the correct file paths and file names
ds = pd.read_excel(os.path.join(path_to_files, metadata_file))
ds["file_path"] = [
    os.path.join(path_to_files, subpath, str(int(session)))
    if pd.notnull(session)
    else None
    for session in ds.session
]
ds["file_name"] = ds.name + ".old"
ds["file_good"] = ds.file_name.notnull()

# Define some settings we need to read your data files properly
kwargs_tiamo = dict(
    file_type="tiamo_de",
    solve_mode="pH_gran",  # can try "pH" or "pH_adjust" or "pH_gran"
)

# First, manually set the acid concentration (`titrant_molinity`) and solve
ds["titrant_molinity"] = titrant_molinity_guess
calk.solve(ds, **kwargs_tiamo)
ds["alkalinity_uncalibrated"] = ds.alkalinity.copy()

# Then, calibrate with CRMs and solve for alkalinity again
calk.calibrate(ds, **kwargs_tiamo)

# Export results to a new spreadsheet
ds.to_excel(
    os.path.join(
        path_to_files,
        metadata_file.replace(".xlsx", " processed.xlsx"),
    ),
    index=False,
)
# Columns of interest
# -------------------
#   alkalinity_uncalibrated
#     TA with the manually set acid concentration (µmol/kg-sw)
#   alkalinity
#     Final TA values, calibrated to CRMs (µmol/kg-sw)
#   alkalinity_offset
#     Difference between measured and certified CRM values (µmol/kg-sw)
#   titrant_molinity_here
#     Best-fitting acid concentration for each CRM (mol/kg)
#   titrant_molinity
#     Average of titrant_molinity_here, used for the final alkalinity values
