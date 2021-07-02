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

# Set MCS parameters and solve with PyCO2SYS
pH = np.append(rng.uniform(size=npts - 1, low=3, high=10), 8.1)
dic = np.append(rng.uniform(size=npts - 1, low=0, high=5000), 0)
results_pyco2 = pyco2.sys(pH, dic, 3, 2, **pyco2_kwargs)
alkalinity_pyco2 = results_pyco2["alkalinity"]

# Solve with Calkulate
results_for_calk = {
    k: v * 1e-6 if k.startswith("total_") or k == "dic" else v
    for k, v in results_pyco2.items()
}
components_calk = calk.simulate.alkalinity_components(
    pH, results_for_calk, results_for_calk
)
alkalinity_calk = calk.simulate.alkalinity(pH, results_for_calk, results_for_calk) * 1e6

# Compare individual components
results_pyco2["H"] = results_pyco2["Hfree"]
results_pyco2["alk_alpha"] = results_pyco2["alkalinity_alpha"]
results_pyco2["alk_beta"] = results_pyco2["alkalinity_beta"]
compare_components = {
    k: np.array([results_pyco2[k], components_calk[k] * 1e6])
    for k in calk.simulate.component_multipliers
    if k != "alkalinity_estimate"
}


def test_components():
    """Does Calkulate's alkalinity components simulator agree with PyCO2SYS?"""
    for v in compare_components.values():
        assert np.all(np.isclose(v[0], v[1], rtol=0, atol=1e-10))


def test_alkalinity_from_pH():
    """Does Calkulate's total alkalinity simulator agree with PyCO2SYS?"""
    assert np.all(np.isclose(alkalinity_calk, alkalinity_pyco2, rtol=0, atol=1e-10))


# test_components()
# test_alkalinity_from_pH()
