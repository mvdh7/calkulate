import calkulate as calk, numpy as np

# Import file
file_name = "tests/data/seawater-CRM-144.dat"
titrant_volume, emf, temperature = calk.io.read_dat(file_name)
titrant_mass = titrant_volume * calk.density.HCl_NaCl_25C_DSC07() * 1e-3
analyte_mass = 0.1  # kg
dic = 2121  # micromol/kg-solution
salinity = 34.1  # practical
nutrients = {
    "total_silicate": 100,
    "total_phosphate": 20,
    "total_ammonia": 12,
    "total_sulfide": 13,
    "total_alpha": 25,
    "total_beta": 25,
}  # all micromol/kg-solution
k_alpha = 1e-5
k_beta = 1e-6


def test_imported_file():
    """Was the titration file imported successfully?"""
    assert isinstance(titrant_volume, np.ndarray)
    assert isinstance(emf, np.ndarray)
    assert isinstance(temperature, np.ndarray)
    assert np.shape(titrant_volume) == np.shape(emf) == np.shape(temperature)


# Get totals and k_constants
totals, totals_pyco2 = calk.interface.get_totals(salinity, dic=dic, **nutrients)
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
alkalinity_certified = 2345  # micromol/kg-solution
solver_args = (titrant_mass, emf, temperature, analyte_mass, totals, k_constants)
res_titrant_molinity = calk.core.calibrate(
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
opt_result = calk.core.solve_emf_complete(
    titrant_molinity,
    titrant_mass,
    emf,
    temperature,
    analyte_mass,
    totals,
    k_constants,
)
alkalinity, emf0 = opt_result["x"]
alkalinity *= 1e6

# Solve for alkalinity with self-calibrated titrant_molinity and user's EMF0 guess
opt_result__emf0_user = calk.core.solve_emf_complete(
    titrant_molinity,
    titrant_mass,
    emf,
    temperature,
    analyte_mass,
    totals,
    k_constants,
    emf0_guess=emf0 + 10,
)
alkalinity__emf0_user, emf0__emf0_user = opt_result["x"]
alkalinity__emf0_user *= 1e6


def test_self_calibration():
    """Does the self-calibrated sample have the certified alkalinity value?"""
    assert np.isclose(alkalinity, alkalinity_certified, rtol=0, atol=1e-6)
    assert 600 < emf0 < 700  # a sensible range for the test file


def test_user_emf0():
    """Do we get the solution with a user-provided first EMF0 guess?"""
    assert np.isclose(alkalinity, alkalinity__emf0_user, rtol=0, atol=1e-6)
    assert np.isclose(emf0, emf0__emf0_user, rtol=0, atol=1e-6)


# Compare with calk.titration functions
ctf_kwargs = dict(
    analyte_mass=analyte_mass,
    dic=dic,
    **nutrients,
    k_alpha=k_alpha,
    k_beta=k_beta,
)
titrant_molinity__ctf, analyte_mass__ctf = calk.titration.calibrate(
    file_name, salinity, alkalinity_certified, **ctf_kwargs
)
(
    alkalinity__ctf,
    emf0__ctf,
    pH__ctf,
    temperature0__ctf,
    analyte_mass__ctf,
    opt_result__ctf,
) = calk.titration.solve(file_name, salinity, titrant_molinity__ctf, **ctf_kwargs)


def test_calibrate_titration():
    """Do the calk.titration functions give consistent results?"""
    assert np.isclose(titrant_molinity, titrant_molinity__ctf, rtol=0, atol=1e-12)
    assert np.isclose(alkalinity, alkalinity__ctf, rtol=0, atol=1e-12)
    assert np.isclose(emf0, emf0__ctf, rtol=0, atol=1e-12)


# test_imported_file()
# test_self_calibration()
# test_user_emf0()
# test_calibrate_titration()
