# %%
import matplotlib.pyplot as plt
import numpy as np
import PyCO2SYS as pyco2
from scipy.optimize import least_squares

from calkulate.convert import emf_to_pH, get_dilution_factor, pH_to_emf
from calkulate.core import (
    _lsqfun_solve_emf_complete,
    calibrate_pH_adjust,
    solve_emf_pH_adjust,
)
from calkulate.read import read_tiamo_de
from calkulate.titration import get_totals_k_constants


volume, pH, temperature = read_tiamo_de("tests/data/ts-tiamo/CRM-210-0033.old")
emf0_offset = -9.830345
emf = pH_to_emf(pH, 0, temperature)
pH = emf_to_pH(emf, emf0_offset, temperature)

salinity = 33.231
dic_certified = 2046.37
alkalinity_certified = 2220.62

analyte_volume = 25
dilfac = get_dilution_factor(volume, analyte_volume)

alkalinity = pyco2.sys(
    par1=pH,
    par1_type=3,
    par2=dic_certified * dilfac * 0,
    par2_type=2,
    temperature=temperature,
    salinity=salinity * dilfac,
    opt_pH_scale=3,
)["alkalinity"]


titrant_molinity = 0.01 * 1e6
# alkalinity_init = alkalinity + volume * 230
alkalinity_init = (
    alkalinity + volume * titrant_molinity / (volume + analyte_volume)
) / dilfac

fig, ax = plt.subplots(dpi=300)
ax.scatter(volume, alkalinity_init)

totals, k_constants = get_totals_k_constants(
    volume,
    temperature,
    analyte_volume,
    salinity,
    dic=0,
)
cal = calibrate_pH_adjust(
    alkalinity_certified,
    volume,
    pH,
    temperature,
    analyte_volume,
    totals,
    k_constants,
    titrant_molinity_guess=0.01 * 1e6,
)

test = solve_emf_pH_adjust(
    titrant_molinity,
    volume * 1e-6,
    pH,
    temperature,
    analyte_volume,
    totals,
    k_constants,
)
print(test["x"][0] * 1e6, test["x"][1])

# Get initial guesses
emf0_guess = 0
# Set which data points to use in the final solver
G = (pH >= 3) & (pH <= 4)
totals_G = {k: v[G] if np.size(v) > 1 else v for k, v in totals.items()}
k_constants_G = {
    k: v[G] if np.size(v) > 1 else v for k, v in k_constants.items()
}
# Solve for alkalinity and EMF0
opt_result = least_squares(
    _lsqfun_solve_emf_complete,
    [alkalinity_certified * 1e-6, emf0_guess],
    args=(
        0.01 * 0.992311 * 1e6,
        volume[G] * 1e-6,
        emf[G],
        temperature[G],
        analyte_volume,
        totals_G,
        k_constants_G,
    ),
    x_scale=[1e-6, 1],
    # **least_squares_kwargs,
)
print(opt_result["x"][0] * 1e6, opt_result["x"][1])
# # Add which data points were used and initial guesses to the output
# opt_result["data_used"] = G
# opt_result["alkalinity_gran"] = alkalinity_guess * 1e6
# opt_result["emf0_gran"] = emf0_guess
