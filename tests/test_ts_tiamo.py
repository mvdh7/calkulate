# %%
from os import listdir

import calkulate as calk


filepath = "tests/data/ts-tiamo/"
filenames = [f for f in listdir(filepath) if f.endswith(".old")]


def test_read_dat():
    for filename in filenames:
        td = calk.read.read_tiamo_de(filepath + filename)
        assert isinstance(td, calk.read.TiamoData)
        assert all(td[0] == td.volume)
        assert all(td[1] == td.pH)
        assert all(td[2] == td.temperature)
        assert td.volume[0] > 0


def test_get_dat_data():
    for filename in filenames:
        td = calk.titration.get_dat_data(
            filepath + filename, read_dat_kwargs={"method": "tiamo_de"}
        )
        assert isinstance(td, calk.titration.DatData)
        assert all(td[0] == td.titrant_mass)
        assert all(td[1] == td.measurement)
        assert all(td[2] == td.temperature)


def test_prepare():
    for filename in filenames:
        pr = calk.titration.prepare(
            filepath + filename,
            35,
            analyte_volume=25,
            kwargs_dat_data={"kwargs_read_dat": {"method": "tiamo_de"}},
        )
        assert isinstance(pr, calk.titration.PrepareResult)
        assert all(pr[0] == pr.titrant_mass)
        assert all(pr[1] == pr.measurement)
        assert all(pr[2] == pr.temperature)
        assert pr[3] == pr.analyte_mass
        assert pr[4] == pr.totals
        assert pr[5] == pr.k_constants
        assert isinstance(pr.totals, dict)
        assert isinstance(pr.k_constants, dict)


# test_read_dat()
# test_get_dat_data()
# test_prepare()
