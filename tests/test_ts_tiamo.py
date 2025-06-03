# %%
import os

import numpy as np

import calkulate as calk


filepath = "tests/data/ts-tiamo/"
filenames = [f for f in os.listdir(filepath) if f.endswith(".old")]


def test_read_dat():
    for filename in filenames:
        dd = calk.read.titrations.read_tiamo_de(
            os.path.join(filepath, filename)
        )
        assert all(dd[0] == dd.titrant_amount)
        assert all(dd[1] == dd.measurement)
        assert all(dd[2] == dd.temperature)
        assert dd.titrant_amount[0] > 0


def test_get_dat_data():
    for filename in filenames:
        dd = calk.read_dat(
            os.path.join(filepath, filename),
            file_type="tiamo_de",
        )
        assert all(dd[0] == dd.titrant_amount)
        assert all(dd[1] == dd.measurement)
        assert all(dd[2] == dd.temperature)


def test_cau():
    for filename in filenames:
        dd = calk.read_dat(
            os.path.join(filepath, filename),
            file_type="tiamo_de",
        )
        cv = calk.convert.amount_units(dd, 33.231, analyte_volume=25)
        assert isinstance(cv, calk.convert.Converted)
        assert all(cv[0] == cv.titrant_mass)
        assert all(cv[1] == cv.measurement)
        assert all(cv[2] == cv.temperature)
        assert cv[3] == cv.analyte_mass


def test_calibrate_solve_emf():
    alkalinity_certified = 2220.62
    salinity = 33.231
    for filename in filenames:
        dd = calk.read_dat(
            os.path.join(filepath, filename),
            file_type="tiamo_de",
        )
        cv = calk.convert.amount_units(dd, salinity, analyte_volume=25)
        emf = calk.convert.pH_to_emf(cv.measurement, 0, cv.temperature)
        totals, k_constants = calk.core.totals_ks(cv)
        cal = calk.core.calibrate_emf(
            alkalinity_certified,
            cv.titrant_mass,
            emf,
            cv.temperature,
            cv.analyte_mass,
            totals,
            k_constants,
            emf0_initial=0,
            titrant_molinity_initial=0.01,
        )
        titrant_molinity = cal["x"][0]
        assert isinstance(titrant_molinity, float)
        assert 0.0095 < titrant_molinity < 0.0100
        sr = calk.core.solve_emf(
            titrant_molinity,
            cv.titrant_mass,
            emf,
            cv.temperature,
            cv.analyte_mass,
            totals,
            k_constants,
            emf0_initial=0,
        )
        assert np.abs(sr.emf0) < 3
        assert np.isclose(sr.alkalinity, alkalinity_certified)


def test_calibrate_solve_pH():
    alkalinity_certified = 2220.62
    salinity = 33.231
    for filename in filenames:
        dd = calk.read_dat(
            os.path.join(filepath, filename),
            file_type="tiamo_de",
        )
        cv = calk.convert.amount_units(dd, salinity, analyte_volume=25)
        totals, k_constants = calk.core.totals_ks(cv, opt_pH_scale=3)
        cal = calk.core.calibrate_pH(
            alkalinity_certified,
            cv.titrant_mass,
            cv.measurement,
            cv.temperature,
            cv.analyte_mass,
            totals,
            k_constants,
            pH_min=3,
            pH_max=4,
            titrant_molinity_initial=0.1,
            titrant_normality=1,
        )
        titrant_molinity = cal["x"][0]
        assert isinstance(titrant_molinity, float)
        assert 0.0095 < titrant_molinity < 0.0100
        spr = calk.core.solve_pH(
            titrant_molinity,
            cv.titrant_mass,
            cv.measurement,
            cv.temperature,
            cv.analyte_mass,
            totals,
            k_constants,
            pH_min=3,
            pH_max=4,
            titrant_normality=1,
        )
        assert np.isclose(spr.alkalinity, alkalinity_certified)


# test_read_dat()
# test_get_dat_data()
# test_cau()
# test_calibrate_solve_emf()
# test_calibrate_solve_pH()
