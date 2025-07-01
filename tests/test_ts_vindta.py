# %%
import numpy as np

import calkulate as calk


# def test_ts_vindta():
dbs = calk.read_dbs(
    "tests/data/ts-vindta/TA_2023.dbs",
    file_path="tests/data/ts-vindta",
    filename_format="{s}-{c} {n} ({d}){b}.dat",
)
dbs.loc[dbs.run_type == "CRM", "total_silicate"] = 3.5
dbs.loc[dbs.run_type == "CRM", "total_phosphate"] = 0.47
dbs.loc[dbs.run_type == "CRM", "dic"] = 2033.86
# ^ should provide at least approx. DIC value for samples too
calk.calibrate(dbs)
# Check values from v23.7.0
assert np.allclose(
    dbs.alkalinity,
    [
        2233.923911,
        2218.556823,
        2218.243184,
    ],
)


# test_ts_vindta()
