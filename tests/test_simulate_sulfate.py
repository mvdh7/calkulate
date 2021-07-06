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
    titrant_mass * 1e3,  # g
    emf,  # mV
    temperature,  # °C
    mode="w",
    measurement_fmt=".4f",
)
file_name_v = "tests/data/test_simulate_sulfate_ml.dat"
titrant_volume = titrant_mass / calk.density.H2SO4_25C_EAIM(titrant_molinity)
calk.write_dat(
    file_name_v,
    titrant_volume * 1e3,  # ml
    emf,  # mV
    temperature,  # °C
    mode="w",
    titrant_amount_fmt=".8f",
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
    titrant_molinity,
    titrant_mass,
    emf,
    temperature,
    analyte_mass,
    totals,
    k_constants,
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
    titrant="H2SO4",
    titrant_amount_unit="g",
)
titrant_molinity_tcal = calk.titration.calibrate(
    file_name,
    salinity,
    alkalinity_core,
    **prepare_kwargs,
)[0]
alkalinity_tcal, emf0_tcal, pH_initial_tcal = calk.titration.solve(
    file_name,
    salinity,
    titrant_molinity_tcal,
    **prepare_kwargs,
)[:3]

# Now volume version
prepare_kwargs_v = dict(
    analyte_mass=analyte_mass,
    analyte_total_sulfate=kwargs_titration["total_sulfate"][0],
    dic=dic,
    molinity_H2SO4=titrant_molinity,
    titrant="H2SO4",
    titrant_amount_unit="ml",
)
titrant_molinity_tcal_v = calk.titration.calibrate(
    file_name_v,
    salinity,
    alkalinity_core,
    **prepare_kwargs_v,
)[0]
alkalinity_tcal_v, emf0_tcal_v, pH_initial_tcal_v = calk.titration.solve(
    file_name_v,
    salinity,
    titrant_molinity_tcal_v,
    **prepare_kwargs_v,
)[:3]

# Import as a Calkulate Dataset, self-calibrate and solve
ds = calk.Dataset({"file_name": [file_name]})
ds["salinity"] = co2sys_core["salinity"]
ds["alkalinity_certified"] = alkalinity_core
ds["analyte_mass"] = analyte_mass
ds["titrant_amount_unit"] = "g"
ds["opt_total_borate"] = 1
ds["opt_k_carbonic"] = 16
ds["dic"] = co2sys_core["dic"]
ds["titrant"] = "H2SO4"
ds.calkulate()

# And again for the volume version
ds_v = calk.Dataset({"file_name": [file_name_v]})
ds_v["salinity"] = co2sys_core["salinity"]
ds_v["alkalinity_certified"] = alkalinity_core
ds_v["analyte_mass"] = analyte_mass
ds_v["titrant_amount_unit"] = "ml"
ds_v["titrant_density"] = np.nan
ds_v["temperature_override"] = np.nan
ds_v["opt_total_borate"] = 1
ds_v["opt_k_carbonic"] = 16
ds_v["dic"] = co2sys_core["dic"]
ds_v["titrant"] = "H2SO4"
ds_v["molinity_H2SO4"] = titrant_molinity
ds_v.calkulate()


def test_calibrate_H2SO4():
    """Do the H2SO4 calibrators find the correct titrant_molinity for a simulated
    titration?
    """
    assert np.isclose(titrant_molinity, titrant_molinity_calibrated, rtol=0, atol=1e-12)
    # These are negligibly worse due to rounding errors when data are saved to file:
    assert np.isclose(titrant_molinity, titrant_molinity_tcal, rtol=0, atol=1e-6)
    assert np.isclose(titrant_molinity, titrant_molinity_tcal_v, rtol=0, atol=1e-6)
    assert np.isclose(titrant_molinity, ds.titrant_molinity_here, rtol=0, atol=1e-6)
    assert np.isclose(titrant_molinity, ds_v.titrant_molinity_here, rtol=0, atol=1e-6)


def test_solve_H2SO4():
    """Do the H2SO4 solvers correctly solve a simulated titration?"""
    assert np.isclose(alkalinity_core, alkalinity_solved, rtol=0, atol=1e-12)
    assert np.isclose(alkalinity_core, alkalinity_cal_solved, rtol=0, atol=1e-12)
    assert np.isclose(alkalinity_core, alkalinity_tcal, rtol=0, atol=1e-10)
    assert np.isclose(alkalinity_core, alkalinity_tcal_v, rtol=0, atol=1e-10)
    assert np.isclose(alkalinity_core, ds.alkalinity, rtol=0, atol=1e-10)
    assert np.isclose(alkalinity_core, ds_v.alkalinity, rtol=0, atol=1e-10)
    assert np.isclose(emf0, emf0_solved, rtol=0, atol=1e-12)
    assert np.isclose(emf0, emf0_cal_solved, rtol=0, atol=1e-12)
    # These are negligibly worse due to rounding errors when data are saved to file:
    assert np.isclose(emf0, emf0_tcal, rtol=0, atol=1e-4)
    assert np.isclose(emf0, emf0_tcal_v, rtol=0, atol=1e-4)
    assert np.isclose(emf0, ds.emf0, rtol=0, atol=1e-4)
    assert np.isclose(emf0, ds_v.emf0, rtol=0, atol=1e-4)
    assert np.isclose(pH_free[0], pH_initial_tcal, rtol=0, atol=1e-6)
    assert np.isclose(pH_free[0], pH_initial_tcal_v, rtol=0, atol=1e-6)
    assert np.isclose(pH_free[0], ds.pH_initial, rtol=0, atol=1e-6)
    assert np.isclose(pH_free[0], ds_v.pH_initial, rtol=0, atol=1e-6)
    assert np.isclose(0, ds.alkalinity_offset, rtol=0, atol=1e-10)
    assert np.isclose(0, ds_v.alkalinity_offset, rtol=0, atol=1e-10)


test_calibrate_H2SO4()
test_solve_H2SO4()
