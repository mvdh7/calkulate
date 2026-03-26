# %%
import os
import shutil

import numpy as np

import calkulate as calk


# Import dbs file
file_path = "tests/data/r2co2/Nico"
dbs = calk.read_dbs("tests/data/r2co2/Nico.dbs", file_path=file_path)

# For each titration file, if there isn't already a backup copy:
# 1. Make a backup copy (with extension .bak).
# 2. Open the .dat, replace first column with correct volumes, and save.
# (because R2-CO2 just saves zeroes in the titrant_amount column!)
tfiles = os.listdir(file_path)
for i, row in dbs.iterrows():
    datfile = os.path.join(row.file_path, row.file_name)
    bakfile = os.path.join(row.file_path, row.file_name)[:-3] + "bak"
    if row.file_name not in tfiles:
        shutil.copyfile(datfile, bakfile)
        dat_data = calk.read_dat(bakfile)
        calk.write_dat(
            datfile,
            # The line below sets the correct titrant_amount values
            np.arange(0, len(dat_data.titrant_amount) * 0.15, 0.15),
            # -----------------------------------------------------
            dat_data.measurement,
            dat_data.temperature,
            mode="w",
            line0=f"Based on '{bakfile}'",
        )

# Now we can proceed with Calkulate as normal
dbs["titrant_molinity"] = 0.1  # guess for now, calibrate later
dbs["salinity"] = 30  # guess for now, measure later
dbs["dic"] = 2350  # guess for now, measure later
dbs["temperature_override"] = 25  # as the files record 0 °C

# Set bad files to ignore
dbs["file_good"] = ~dbs.bottle.isin(["junk-250916-01"])
dbs.solve()

# TODO Workaround for density storing bug in v23.7.0
dbs["analyte_mass"] = (
    1e-3
    * dbs.analyte_volume
    / calk.density.seawater_1atm_MP81(
        temperature=dbs.temperature_override, salinity=dbs.salinity
    )
)

# Extract one titration to plot
tt = dbs.to_Titration(1)
tt.plot_pH()
tt.plot_components()
tt.plot_alkalinity()
