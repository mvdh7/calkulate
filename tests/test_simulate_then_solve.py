# %%
import copy

import numpy as np
import pandas as pd
import PyCO2SYS as pyco2

import calkulate as calk


# Function inputs
dic = np.array([2000.0])
pH_free = np.array([8.1])
temperature = 25.0
salinity = 35.0
# For the titration
analyte_mass = 0.2  # kg
titrant_molinity = 0.3  # mol/kg
titrant_mass_start = 0
titrant_mass_step = 0.05e-3
titrant_mass_stop = 2.51e-3
titrant_mass = np.arange(
    titrant_mass_start, titrant_mass_stop, titrant_mass_step
)  # kg
emf0 = 300  # mV
# ===============

# Assemble inputs into dicts
if np.isscalar(temperature):
    temperature = np.array([temperature])
if np.isscalar(salinity):
    salinity = np.array([salinity])
pyco2sys_kwargs = dict(
    opt_k_carbonic=10,
    opt_total_borate=1,
)
kwargs_core = dict(
    temperature=temperature,
    salinity=salinity,
    opt_pH_scale=3,
    **pyco2sys_kwargs,
)

# Calculate total alkalinity
co2sys_core = pyco2.sys(dic, pH_free, 2, 3, **kwargs_core)
alkalinity_core = co2sys_core["alkalinity"]

# Simulate titration(s)
dilution_factor = analyte_mass / (analyte_mass + titrant_mass)
alkalinity_titration = (
    1e6
    * (
        analyte_mass * co2sys_core["alkalinity"] * 1e-6
        - titrant_mass * titrant_molinity
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
    "total_sulfate",
]:
    kwargs_titration[k] = kwargs_titration[k] * dilution_factor
co2sys_titrations = pyco2.sys(
    alkalinity_titration, dic_titration, 1, 2, **kwargs_titration
)
pH_titrations = co2sys_titrations["pH_free"]

# Export .dat file(s) for Calkulate
emf = calk.convert.pH_to_emf(
    pH_titrations, emf0, kwargs_titration["temperature"]
)
file_name = "tests/data/test_simulate_then_solve.dat"
calk.write_dat(
    file_name,
    titrant_mass * 1000,  # g
    emf,  # mV
    co2sys_titrations["temperature"],  # °C
    mode="w",
    measurement_fmt=".4f",
)

# Test simulate._titration function
(
    titrant_mass__simfunc,
    emf__simfunc,
    temperature__simfunc,
    analyte_mass__simfunc,
    totals__simfunc,
    k_constants__simfunc,
) = calk.simulate._titration(
    alkalinity_core,
    dic=dic,
    emf0=emf0,
    salinity=salinity,
    analyte_mass=analyte_mass,
    temperature=temperature,
    titrant_molinity=titrant_molinity,
    titrant_mass_start=0,
    titrant_mass_step=0.05e-3,
    titrant_mass_stop=2.51e-3,
    **pyco2sys_kwargs,
)


def test_manual_solve():
    # Solve manually
    cv = calk.convert.Converted(
        titrant_mass, emf, temperature, analyte_mass, salinity
    )
    totals, k_constants = calk.core.totals_ks(
        cv,
        dic=dic,
        **pyco2sys_kwargs,
        dilute_totals_for_ks=True,
    )
    sr = calk.core.solve_emf(
        titrant_molinity,
        titrant_mass,
        emf,
        np.full_like(emf, temperature),
        analyte_mass,
        totals,
        k_constants,
        alkalinity_init=None,
        double=True,
        emf0_init=None,
        pH_min=0,
        pH_max=14,
        titrant_normality=1,
    )
    for k, v in totals.items():
        if k in totals__simfunc:
            assert np.allclose(totals__simfunc[k], v), k
    for k, v in k_constants.items():
        if k in k_constants__simfunc:
            assert np.allclose(k_constants__simfunc[k], v), k
    assert np.isclose(alkalinity_core[0], sr.alkalinity)
    # NOTE The problem was that we've switched to not allowing equilibrium
    # constants to change through the titration due to the pH scale correction
    # changing as sulfate gets diluted!  But only in the solver, and not in the
    # simulations here!  Now fixed with the kwarg dilute_totals_for_ks=True,
    # but that is NOT default behaviour from v23.7


# Import as a Calkulate Dataset
ds = pd.DataFrame({"file_name": [file_name]})
ds["salinity"] = co2sys_core["salinity"]
ds["analyte_mass"] = analyte_mass
ds["titrant_molinity"] = titrant_molinity
ds["titrant_amount_unit"] = "g"
ds["opt_total_borate"] = 1
ds["opt_k_carbonic"] = 10
ds["dic"] = co2sys_core["dic"]
ds = calk.Dataset(ds)
ds["dilute_totals_for_ks"] = True
ds.solve()
co2sys_core["alkalinity_titration"] = alkalinity_solved = (
    ds.alkalinity.to_numpy()[0]
)
co2sys_core["emf0"] = ds.emf0.to_numpy()
co2sys_core["alkalinity_offset"] = (
    co2sys_core["alkalinity_titration"] - alkalinity_core
)


def test_simulate_then_solve():
    assert np.isclose(co2sys_core["emf0"], emf0, rtol=0, atol=1e-4)
    assert np.isclose(co2sys_core["alkalinity_offset"], 0, rtol=0, atol=1e-3)


def test_simulate_titration_function():
    assert np.allclose(titrant_mass, titrant_mass__simfunc)
    assert np.allclose(emf, emf__simfunc)
    assert np.allclose(temperature, temperature__simfunc)


# test_manual_solve()
# test_simulate_then_solve()
# test_simulate_titration_function()
