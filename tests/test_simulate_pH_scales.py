import numpy as np
import PyCO2SYS as pyco2
import calkulate as calk

# Set conditions
salinity = 35
temperature = 25
dic = 2100
alkalinity_initial = -500

# Create empty dicts
pH = {}
k_constants = {}
components = {}
alkalinity = {}
all_H = {}

# Get pH on all scales
res = pyco2.sys(
    par1=alkalinity_initial,
    par1_type=1,
    par2=dic,
    par2_type=2,
    temperature=temperature,
    salinity=salinity,
)
pH[1] = res["pH_total"]
pH[2] = res["pH_sws"]
pH[3] = res["pH_free"]

# Get totals (pH-scale independent)
totals, totals_pyco2 = calk.interface.get_totals(
    salinity,
    dic=dic,
)

for i in pH.keys():
    # Get k_constants on the pH scale
    k_constants[i] = calk.interface.get_k_constants(
        totals_pyco2,
        temperature,
        opt_pH_scale=i,
    )

    # Simulate solution_components
    components[i] = calk.simulate.alkalinity_components(
        pH[i], totals, k_constants[i], opt_pH_scale=i
    )
    alkalinity[i] = (
        calk.simulate.alkalinity(pH[i], totals, k_constants[i], opt_pH_scale=i) * 1e6
    )
    all_H[i] = components[i]["H"]
    if "HSO4" in components[i]:
        all_H[i] += components[i]["HSO4"]
    if "HF" in components[i]:
        all_H[i] += components[i]["HF"]
    all_H[i] *= 1e6


def test_nonH_components():
    """Are all non-H components identical regardless of the pH scale?"""
    nonH_components = ["BOH4", "CO2", "HCO3", "CO3", "OH"]
    tol = dict(atol=0, rtol=1e-12)
    for nonH in nonH_components:
        for i, j in ((1, 2), (1, 3), (2, 3)):
            assert np.isclose(components[i][nonH], components[j][nonH], **tol)


def test_H_components():
    """Are the H components consistent across the different pH scales?"""
    tol = dict(atol=0, rtol=0.005)
    for i, j in ((1, 2), (1, 3), (2, 3)):
        assert np.isclose(all_H[i], all_H[j], **tol)


def test_alkalinity():
    """Is the alkalinity consistent across the different pH scales?"""
    tol = dict(atol=0, rtol=0.005)
    for i in [1, 2, 3]:
        assert np.isclose(alkalinity[i], alkalinity_initial, **tol)


# test_nonH_components()
# test_H_components()
# test_alkalinity()
