import copy
import numpy as np
import PyCO2SYS as pyco2, calkulate as calk


# Function inputs
dic = np.array([2000.0])
pH_free = np.array([8.1])
salinity = 35.0
# For the titration
analyte_mass = 0.2  # kg
titrant_molinity = 0.15  # mol/kg
titrant_mass = np.arange(0, 2.51, 0.05) * 1e-3  # kg
temperature = np.full_like(titrant_mass, 25.0)
titrant_alkalinity_factor = 2  # for H2SO4
emf0 = 300  # mV
# ===============

# Assemble inputs into dicts
if np.isscalar(temperature):
    temperature = np.array([temperature])
if np.isscalar(salinity):
    salinity = np.array([salinity])
kwargs_core = dict(
    temperature=temperature[0],
    salinity=salinity,
    opt_pH_scale=3,
    opt_k_carbonic=calk.default.opt_k_carbonic,
    opt_total_borate=calk.default.opt_total_borate,
)

# Calculate total alkalinity
co2sys_core = pyco2.sys(dic, pH_free, 2, 3, **kwargs_core)
alkalinity_core = co2sys_core["alkalinity"]

# Dilute things with the titrant
dilution_factor = analyte_mass / (analyte_mass + titrant_mass)
alkalinity_titration = (
    1e6
    * (
        analyte_mass * co2sys_core["alkalinity"] * 1e-6
        - titrant_alkalinity_factor * titrant_mass * titrant_molinity
    )
    / (analyte_mass + titrant_mass)
)
dic_titration = dic * dilution_factor
kwargs_titration = copy.deepcopy(kwargs_core)
kwargs_titration.update(
    total_borate=co2sys_core["total_borate"],
    total_fluoride=co2sys_core["total_fluoride"],
    total_sulfate=co2sys_core["total_sulfate"],
)
for k in [
    "total_borate",
    "total_fluoride",
]:
    kwargs_titration[k] = kwargs_titration[k] * dilution_factor

# Sulfate gets added by the H2SO4 titrant, not diluted
kwargs_titration["total_sulfate"] = (
    kwargs_titration["total_sulfate"] * analyte_mass
    + titrant_molinity * titrant_mass * 1e6
) / (analyte_mass + titrant_mass)

# Simulate titration with PyCO2SYS
co2sys_titrations = pyco2.sys(
    alkalinity_titration, dic_titration, 1, 2, **kwargs_titration
)
pH_titrations = co2sys_titrations["pH_free"]

# Export .dat file(s) for Calkulate
emf = calk.convert.pH_to_emf(pH_titrations, emf0, kwargs_titration["temperature"])
file_name = "tests/data/test_simulate_sulfate.dat"
calk.write_dat(
    file_name,
    titrant_mass * 1000,  # g
    emf,  # mV
    temperature,  # °C
    mode="w",
    measurement_fmt=".4f",
)

# Get totals and k_constants
totals, totals_pyco2 = calk.interface.get_totals(salinity, dic=dic)
totals = calk.convert.dilute_totals_H2SO4(
    totals, titrant_molinity, titrant_mass, analyte_mass
)
totals_pyco2 = calk.convert.totals_to_pyco2(totals, salinity)
k_constants = calk.interface.get_k_constants(totals_pyco2, temperature)

# Calibrate!
opt_result_calibrate = calk.core.calibrate_H2SO4(
    alkalinity_core,
    titrant_mass,
    emf,
    temperature,
    analyte_mass,
    kwargs_titration["total_sulfate"][0] * 1e-6,
    salinity,
    totals,
    k_constants,
)
titrant_molinity_calibrated = opt_result_calibrate["x"][0]

# Solve!
opt_result_solve = calk.core.solve_emf_complete_H2SO4(
    titrant_molinity, titrant_mass, emf, temperature, analyte_mass, totals, k_constants,
)
alkalinity_solved, emf0_solved = opt_result_solve["x"]
alkalinity_solved *= 1e6

# And again with the calibrated titrant_molinity
opt_result_cal_solve = calk.core.solve_emf_complete_H2SO4(
    titrant_molinity_calibrated,
    titrant_mass,
    emf,
    temperature,
    analyte_mass,
    totals,
    k_constants,
)
alkalinity_cal_solved, emf0_cal_solved = opt_result_cal_solve["x"]
alkalinity_cal_solved *= 1e6

# Try titration calibrate/solve
prepare_kwargs = dict(
    analyte_mass=analyte_mass,
    analyte_total_sulfate=kwargs_titration["total_sulfate"][0],
    dic=dic,
    molinity_H2SO4=titrant_molinity,
    titrant="H2SO4",
    titrant_amount_unit="g",
)
titrant_molinity_tcal = calk.titration.calibrate(
    file_name, salinity, alkalinity_core, **prepare_kwargs,
)[0]
alkalinity_tcal, emf0_tcal, pH_initial_tcal = calk.titration.solve(
    file_name, salinity, titrant_molinity_tcal, **prepare_kwargs,
)[:3]


def test_calibrate_H2SO4():
    """Do the H2SO4 calibrators find the correct titrant_molinity for a simulated
    titration?
    """
    assert np.isclose(titrant_molinity, titrant_molinity_calibrated, rtol=0, atol=1e-12)
    assert np.isclose(titrant_molinity, titrant_molinity_tcal, rtol=0, atol=1e-6)
    # ^ this one is negligibly worse due to rounding errors when data are saved to file


def test_solve_H2SO4():
    """Do the H2SO4 solvers correctly solve a simulated titration?"""
    assert np.isclose(alkalinity_core, alkalinity_solved, rtol=0, atol=1e-12)
    assert np.isclose(emf0, emf0_solved, rtol=0, atol=1e-12)
    assert np.isclose(alkalinity_core, alkalinity_cal_solved, rtol=0, atol=1e-12)
    assert np.isclose(emf0, emf0_cal_solved, rtol=0, atol=1e-12)
    assert np.isclose(alkalinity_core, alkalinity_tcal, rtol=0, atol=1e-12)
    assert np.isclose(emf0, emf0_tcal, rtol=0, atol=1e-4)
    # ^ this one is negligibly worse due to rounding errors when data are saved to file
    assert np.isclose(pH_free[0], pH_initial_tcal, rtol=0, atol=1e-6)


test_calibrate_H2SO4()
test_solve_H2SO4()

# # Import as a Calkulate Dataset
# ds = calk.Dataset({"file_name": [file_name]})
# ds["salinity"] = co2sys_core["salinity"]
# ds["analyte_mass"] = analyte_mass
# ds["titrant_molinity"] = titrant_molinity
# ds["titrant_amount_unit"] = "g"
# ds["opt_total_borate"] = 1
# ds["opt_k_carbonic"] = 16
# ds["dic"] = co2sys_core["dic"]
# ds["titrant"] = "H2SO4"
# ds.solve()
# co2sys_core["alkalinity_titration"] = alkalinity_solved = ds.alkalinity.to_numpy()[0]
# co2sys_core["emf0"] = ds.emf0.to_numpy()
# co2sys_core["alkalinity_offset"] = co2sys_core["alkalinity_titration"] - alkalinity_core


# def test_simulate_then_solve():
#     assert np.isclose(co2sys_core["emf0"], emf0, rtol=0, atol=1e-4)
#     assert np.isclose(co2sys_core["alkalinity_offset"], 0, rtol=0, atol=1e-3)


# # test_simulate_then_solve()
