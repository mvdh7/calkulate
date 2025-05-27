# %%
from os import listdir

import numpy as np

import calkulate as calk


filepath = "tests/data/ts-tiamo/"
filenames = [f for f in listdir(filepath) if f.endswith(".old")]


def test_read_dat():
    for filename in filenames:
        dd = calk.read.titrations.read_tiamo_de(filepath + filename)
        assert all(dd[0] == dd.titrant_amount)
        assert all(dd[1] == dd.measurement)
        assert all(dd[2] == dd.temperature)
        assert dd.titrant_amount[0] > 0


def test_get_dat_data():
    for filename in filenames:
        dd = calk.read_dat(filepath + filename, file_type="tiamo_de")
        assert all(dd[0] == dd.titrant_amount)
        assert all(dd[1] == dd.measurement)
        assert all(dd[2] == dd.temperature)


def test_cau():
    for filename in filenames:
        dd = calk.read_dat(filepath + filename, file_type="tiamo_de")
        cv = calk.titration.convert_amount_units(dd, 33.231, analyte_volume=25)
        assert isinstance(cv, calk.titration.Converted)
        assert all(cv[0] == cv.titrant_mass)
        assert all(cv[1] == cv.measurement)
        assert all(cv[2] == cv.temperature)
        assert cv[3] == cv.analyte_mass


def test_tt_calibrate_solve():
    alkalinity_certified = 2220.62
    for filename in filenames:
        dd = calk.read_dat(filepath + filename, file_type="tiamo_de")
        cv = calk.titration.convert_amount_units(dd, 33.231, analyte_volume=25)
        totals, k_constants = calk.titration.get_totals_k_constants(cv, 33.231)
        titrant_molinity = calk.titration.calibrate(
            alkalinity_certified,
            cv,
            totals,
            k_constants,
            measurement_type="pH",
            titrant_molinity_guess=0.01,
        )
        assert isinstance(titrant_molinity, float)
        assert 0.0095 < titrant_molinity < 0.0100
        sr = calk.titration.solve(
            titrant_molinity,
            cv,
            totals,
            k_constants,
            measurement_type="pH",
        )
        assert np.isclose(sr.alkalinity, alkalinity_certified)
        assert np.abs(sr.emf0) < 3


# test_read_dat()
# test_get_dat_data()
# test_cau()
# test_tt_calibrate_solve()
