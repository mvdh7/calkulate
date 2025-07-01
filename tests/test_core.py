# %%
import numpy as np

import calkulate as calk


def test_imported_file():
    """Was the titration file imported successfully?"""
    # Import file
    file_name = "tests/data/seawater-CRM-144.dat"
    titrant_volume, emf, temperature = calk.read_dat(file_name)
    assert isinstance(titrant_volume, np.ndarray)
    assert isinstance(emf, np.ndarray)
    assert isinstance(temperature, np.ndarray)
    assert np.shape(titrant_volume) == np.shape(emf) == np.shape(temperature)


def test_self_calibration():
    """Does the self-calibrated sample have the certified alkalinity value?"""
    # NOTE that the code here is old and no longer the most efficient way to
    # get to the input arguments required for the solve and calibrate functions
    # Import file
    file_name = "tests/data/seawater-CRM-144.dat"
    titrant_volume, emf, temperature = calk.read_dat(file_name)
    titrant_mass = titrant_volume * calk.density.HCl_NaCl_25C_DSC07() * 1e-3
    analyte_mass = 0.1  # kg
    dic = 2121  # µmol/kg-sol
    salinity = 34.1  # practical
    nutrients = {
        "total_silicate": 100,
        "total_phosphate": 20,
        "total_ammonia": 12,
        "total_sulfide": 13,
        "total_alpha": 25,
        "total_beta": 25,
    }  # all µmol/kg-sol
    k_alpha = 1e-5
    k_beta = 1e-6

    # Get totals and k_constants
    totals, totals_pyco2 = calk.interface.get_totals(
        salinity, dic=dic, **nutrients
    )
    totals = calk.convert.dilute_totals(totals, titrant_mass, analyte_mass)
    totals_pyco2 = calk.convert.dilute_totals_pyco2(
        totals_pyco2, titrant_mass, analyte_mass
    )
    k_constants = calk.interface.get_k_constants(
        totals_pyco2,
        temperature,
        k_alpha=k_alpha,
        k_beta=k_beta,
    )

    # Solve for titrant_molinity
    alkalinity_certified = 2345  # µmol/kg-sol
    res_titrant_molinity = calk.core.calibrate_emf(
        alkalinity_certified,
        titrant_mass,
        emf,
        temperature,
        analyte_mass,
        totals,
        k_constants,
    )
    titrant_molinity = res_titrant_molinity["x"][0]

    # Solve for alkalinity with self-calibrated titrant_molinity
    sr = calk.core.solve_emf(
        titrant_molinity,
        titrant_mass,
        emf,
        temperature,
        analyte_mass,
        totals,
        k_constants,
    )
    assert np.isclose(sr.alkalinity, alkalinity_certified, rtol=0, atol=1e-6)
    assert 600 < sr.emf0 < 700  # a sensible range for the test file


# test_imported_file()
# test_self_calibration()
