# %%
import numpy as np
import pandas as pd

import calkulate as calk


salinity = 29
analyte_volume = 95  # ml
titrant_molinity = 0.09477679299470473  # µmol/kg

df = {
    "index": np.array(["1a", "1b", "2a", "2b"]),
    "file_name": np.array(
        [
            "12-6  25  (0)STANDARD-1A.dat",
            "4-6  25  (0)STANDARD-1B.dat",
            "4-6  25  (0)STANDARD-2A.dat",
            "4-6  25  (0)STANDARD-2B.dat",
        ]
    ),
    "salinity": salinity,
    "file_path": "tests/data/seao2",
    "analyte_volume": analyte_volume,
    "titrant_molinity": titrant_molinity,
    "dic": [2067.3, 2065.2, 2067.5, 2073.1],
    "opt_k_fluoride": 2,
    # # These aren't strictly needed because they're the defaults, but here
    # # you can see which options were used (see PyCO2SYS docs for details):
    # "opt_k_carbonic": 10,
    # "opt_total_borate": 1,
    # "opt_k_bisulfate": 1,
}
df = pd.DataFrame(df).set_index("index")
calk.solve(
    df,
    # These kwargs revert changes of behaviour to pre-v23.7:
    dilute_totals_for_ks=True,
    double=False,
)

# print(df.alkalinity)
# print(df.alkalinity.mean())
# print(df.alkalinity.std())


def test_seao2_opt1():
    assert (
        np.round(df.alkalinity, 1) == [2132.4, 2132.1, 2133.4, 2132.2]
    ).all()
    assert np.round(df.alkalinity.mean(), 4) == 2132.5225


# test_seao2_opt1()
