import calkulate as calk, PyCO2SYS as pyco2, numpy as np

# Set up random number generator for repeatability
rng = np.random.default_rng(7)

# Define conditions
npts = 10000
pyco2_kwargs = {
    "salinity": rng.uniform(size=npts, low=1, high=50),
    "temperature": rng.uniform(size=npts, low=0, high=50),
    "pressure": rng.uniform(size=npts, low=0, high=10000),
    "total_silicate": rng.uniform(size=npts, low=0, high=100),
    "total_phosphate": rng.uniform(size=npts, low=0, high=50),
    "total_ammonia": rng.uniform(size=npts, low=0, high=50),
    "total_sulfide": rng.uniform(size=npts, low=0, high=50),
    "opt_k_carbonic": rng.integers(size=npts, low=1, high=16, endpoint=True),
    "opt_total_borate": rng.integers(size=npts, low=1, high=2, endpoint=True),
    "opt_pH_scale": 3,
    "opt_k_bisulfate": rng.integers(size=npts, low=1, high=2, endpoint=True),
    "opt_k_fluoride": rng.integers(size=npts, low=1, high=2, endpoint=True),
    "opt_gas_constant": 3,
    "buffers_mode": "none",
    "total_alpha": rng.uniform(size=npts, low=0, high=100),
    "k_alpha": rng.uniform(size=npts, low=2, high=12),
    "total_beta": rng.uniform(size=npts, low=0, high=100),
    "k_beta": rng.uniform(size=npts, low=2, high=12),
}
# Remove Peng options (bad alkalinity equation!)
opt_kc = pyco2_kwargs["opt_k_carbonic"]
pyco2_kwargs["opt_k_carbonic"][opt_kc == 7] = 16

# Evaluate PyCO2SYS intermediates
totals = pyco2.salts.assemble(
    pyco2_kwargs["salinity"],
    pyco2_kwargs["total_silicate"],
    pyco2_kwargs["total_phosphate"],
    pyco2_kwargs["total_ammonia"],
    pyco2_kwargs["total_sulfide"],
    pyco2_kwargs["opt_k_carbonic"],
    pyco2_kwargs["opt_total_borate"],
)
k_constants = pyco2.equilibria.assemble(
    pyco2_kwargs["temperature"],
    pyco2_kwargs["pressure"],
    totals,
    pyco2_kwargs["opt_pH_scale"],
    pyco2_kwargs["opt_k_carbonic"],
    pyco2_kwargs["opt_k_bisulfate"],
    pyco2_kwargs["opt_k_fluoride"],
    pyco2_kwargs["opt_gas_constant"],
)

# Put everything into a DataFrame for Calkulate
totals_calk = {k: v * 1e6 for k, v in totals.items() if k != "Sal"}
totals_calk.update({k: pyco2_kwargs[k] for k in ["total_alpha", "total_beta"]})
k_constants.update({k: pyco2_kwargs[k] for k in ["k_alpha", "k_beta"]})
k_constants.pop("alpha")
k_constants.pop("beta")
titration = calk.Titration({**totals_calk, **k_constants}).rename(
    mapper=calk.convert.pyco2_to_calk, axis=1
)

# Solve the CO2 system both ways (with Calkulate and with PyCO2SYS)
pH = np.append(rng.uniform(size=npts - 1, low=3, high=10), 8.1)
dic = titration["dic"] = np.append(rng.uniform(size=npts - 1, low=0, high=5000), 0)
alkalinity_calk = calk.simulate.alkalinity(pH, titration) * 1e6
alkalinity_pyco2 = pyco2.CO2SYS_nd(pH, dic, 3, 2, **pyco2_kwargs)["alkalinity"]


def test_alkalinity_from_pH():
    """Does Calkulate's alkalinity simulator agree with PyCO2SYS?"""
    assert np.all(np.isclose(alkalinity_calk, alkalinity_pyco2, rtol=0, atol=1e-10))


test_alkalinity_from_pH()
