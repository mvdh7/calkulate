import copy
import numpy as np, pandas as pd
import PyCO2SYS as pyco2, calkulate as calk


# Function inputs
dic = np.array([2000.0])
pH_free = np.array([8.1])
temperature = 25.0
salinity = 35.0
# For the titration
sample_mass = 0.2  # kg
titrant_molinity = 0.3  # mol/kg
titrant_mass = np.arange(0, 2.51, 0.05) * 1e-3  # kg
emf0 = 300  # mV
# ===============

# Assemble inputs into dicts
if np.isscalar(temperature):
    temperature = np.array([temperature])
if np.isscalar(salinity):
    salinity = np.array([salinity])
kwargs_core = dict(
    temperature=temperature,
    salinity=salinity,
    opt_pH_scale=3,
    opt_k_carbonic=16,
    opt_total_borate=1,
)

# Calculate total alkalinity
co2sys_core = pyco2.sys(dic, pH_free, 2, 3, **kwargs_core)
alkalinity_core = co2sys_core["alkalinity"]

# Simulate titration(s)
dilution_factor = sample_mass / (sample_mass + titrant_mass)
alkalinity_titration = (
    1e6
    * (sample_mass * alkalinity_core * 1e-6 - titrant_mass * titrant_molinity)
    / (sample_mass + titrant_mass)
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
    "total_sulfate",
]:
    kwargs_titration[k] = kwargs_titration[k] * dilution_factor
co2sys_titrations = pyco2.sys(
    alkalinity_titration, dic_titration, 1, 2, **kwargs_titration
)
pH_titrations = co2sys_titrations["pH_free"]

# Export .dat file(s) for Calkulate
file_path = "tests/data/"
file_name = "test_pH_scales_simulated_titration.dat"
emf = calk.convert.pH_to_emf(pH_titrations, emf0, kwargs_titration["temperature"])
calk.write_dat(
    file_path + file_name,
    titrant_mass * 1000,  # g
    emf,  # mV
    co2sys_titrations["temperature"],  # °C
    mode="w",
    measurement_fmt=".8f",
)

# Import as a Calkulate Dataset
ds = pd.DataFrame({"file_name": [file_name, file_name, file_name]})
ds["file_path"] = file_path
ds["analyte_mass"] = sample_mass
ds["titrant_molinity"] = titrant_molinity
ds["titrant_amount_unit"] = "g"
ds["opt_total_borate"] = 1
ds["opt_k_carbonic"] = 16
ds["alkalinity_true"] = alkalinity_core[0]
for v in ["salinity", "dic", "pH_free", "pH_total", "pH_sws"]:
    try:
        ds[v] = co2sys_core[v][0]
    except:
        ds[v] = co2sys_core[v]
ds["opt_pH_scale"] = [1, 2, 3]

# Solve and check
calk.dataset.solve(ds)
