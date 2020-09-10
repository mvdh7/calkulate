import calkulate as calk, PyCO2SYS as pyco2, numpy as np

# Define conditions
pyco2_kwargs = {
    "salinity": 35,
    "temperature": 25,
    "pressure": 0,
    "total_silicate": 10,
    "total_phosphate": 1,
    "total_ammonia": 2,
    "total_sulfide": 1.5,
    "opt_k_carbonic": 16,
    "opt_total_borate": 1,
    "opt_pH_scale": 3,
    "opt_k_bisulfate": 1,
    "opt_k_fluoride": 2,
    "opt_gas_constant": 3,
}

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


def test_alkalinity_from_pH():
    npts = 1000
    pH = np.append(np.random.normal(size=npts, loc=8.1, scale=3), 8.1)
    dic = np.append(np.random.normal(size=npts, loc=2000, scale=400), 0)
    dic[dic < 0] = 0
    alkalinity_calk = (
        calk.simulate.alkalinity(pH, totals, k_constants, dic=dic * 1e-6) * 1e6
    )
    # Switch to fixed TA equation in PyCO2SYS before comparing
    pyco2.solve.get.TAfromTCpH = pyco2.solve.get.TAfromTCpH_fixed
    alkalinity_pyco2 = pyco2.CO2SYS_nd(pH, dic, 3, 2, **pyco2_kwargs)["alkalinity"]
    assert np.all(np.isclose(alkalinity_calk, alkalinity_pyco2, rtol=0, atol=1e-3))


test_alkalinity_from_pH()
