# %%
from os import listdir

import calkulate as calk


def test_read_dat():
    filepath = "tests/data/ts-tiamo/"
    filenames = [f for f in listdir(filepath) if f.endswith(".old")]
    for filename in filenames:
        td = calk.read.read_tiamo_de(filepath + filename)
        assert isinstance(td, calk.read.TiamoData)
        assert all(td[0] == td.volume)
        assert all(td[1] == td.pH)
        assert all(td[2] == td.temperature)
        assert td.volume[0] > 0


# test_read_dat()
