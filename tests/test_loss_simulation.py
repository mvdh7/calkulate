import numpy as np
import PyCO2SYS as pyco2
from calkulate import convert, default

alkalinity = 2300
analyte_mass = 0.1
dic = 2000
emf0 = 600
salinity = 35
temperature = 25
titrant_mass_start = 0
titrant_mass_step = 0.15e-3
titrant_mass_stop = 4.2e-3
titrant_molinity = 0.1
fCO2_air = default.fCO2_air
k_dic_loss = 2
pyco2sys_kwargs = {}

titrant_mass = np.arange(titrant_mass_start, titrant_mass_stop, titrant_mass_step)
if np.isscalar(temperature):
    temperature = np.full_like(titrant_mass, temperature)
# Ensure we use Calkulate's default PyCO2SYS options if they're not provided
# within pyco2sys_kwargs
for opt in [
    "opt_gas_constant",
    "opt_k_bisulfate",
    "opt_k_carbonic",
    "opt_k_fluoride",
    "opt_total_borate",
]:
    if opt not in pyco2sys_kwargs:
        pyco2sys_kwargs[opt] = getattr(default, opt)
# Set up dict of start-condition kwargs and then get totals from PyCO2SYS
kwargs_start = dict(
    salinity=salinity,
    temperature=temperature,
    opt_pH_scale=3,
    **pyco2sys_kwargs,
)
co2sys_start = pyco2.sys(**kwargs_start)
# Calculate dilution through the titration, keeping units in µmol/kg
dilution_factor = analyte_mass / (analyte_mass + titrant_mass)
alkalinity_titration = (
    1e6
    * (analyte_mass * alkalinity * 1e-6 - titrant_mass * titrant_molinity)
    / (analyte_mass + titrant_mass)
)
kwargs_titration = kwargs_start.copy()
for k in [
    "total_alpha",
    "total_ammonia",
    "total_beta",
    "total_borate",
    "total_calcium",
    "total_fluoride",
    "total_phosphate",
    "total_silicate",
    "total_sulfate",
    "total_sulfide",
]:
    kwargs_titration[k] = co2sys_start[k] * dilution_factor
# Model DIC loss
if k_dic_loss is not None:
    dic_titration = np.full_like(titrant_mass, np.nan)
    dic_titration[0] = dic
    for i in range(len(dic_titration) - 1):
        kwargs_i = {
            k: v if np.isscalar(v) else v[i] for k, v in kwargs_titration.items()
        }
        fCO2 = pyco2.sys(
            par1=alkalinity_titration[i],
            par2=dic_titration[i],
            par1_type=1,
            par2_type=2,
            **kwargs_i,
        )["fCO2"]
        delta_fCO2 = fCO2 - fCO2_air
        dic_loss = k_dic_loss * delta_fCO2 * titrant_mass_step
        dic_titration[i + 1] = (
            dic_titration[i]
            * convert.get_dilution_factor(
                titrant_mass_step, analyte_mass + titrant_mass[i]
            )
            - dic_loss
        )
else:
    dic_titration = dic * dilution_factor
# Simulate the titration pH and convert it to EMF
co2sys_titration = pyco2.sys(
    alkalinity_titration, dic_titration, 1, 2, **kwargs_titration
)
pH_titration = co2sys_titration["pH_free"]
emf = convert.pH_to_emf(pH_titration, emf0, temperature)
# Get totals (in mol/kg) and k_constants dicts for other Calkulate functions
totals = {k: v * 1e-6 for k, v in co2sys_titration.items() if k.startswith("total_")}
totals["dic"] = co2sys_titration["dic"] * 1e-6
k_constants = {k: v for k, v in co2sys_titration.items() if k.startswith("k_")}
