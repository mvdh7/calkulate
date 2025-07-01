# %%
import matplotlib.pyplot as plt

import calkulate as calk


file_name = "tests/data/seawater-CRM-144.dat"

dd = calk.read_dat(file_name)
titrant_mass = dd.titrant_amount * calk.density.HCl_NaCl_25C_DSC07() * 1e-3
analyte_mass = 0.1  # kg
titrant_molinity = 0.1
titrant_normality = 1
titrant_total_sulfate = 0
# emf0_init = 657.553
emf0_init = None
alkalinity_certified = 2238.60
salinity = 33.571
dic = 2031.53

cv = calk.convert.amount_units(dd, salinity, analyte_mass=analyte_mass)
totals, k_constants = calk.core.totals_ks(cv, dic=dic)

kwargs_files = dict(
    analyte_mass=analyte_mass,
    dic=dic,
    titrant_normality=titrant_normality,
    titrant_total_sulfate=titrant_total_sulfate,
)
fcal = calk.files.calibrate(
    file_name,
    alkalinity_certified,
    salinity,
    **kwargs_files,
)
fsr = calk.files.solve(
    file_name,
    fcal["x"][0],
    salinity,
    **kwargs_files,
)

# Guess alkalinity and EMF0
ggr = calk.core.gran_guesses(
    titrant_mass,
    dd.measurement,
    dd.temperature,
    analyte_mass,
    titrant_molinity,
)
cal = calk.core.calibrate_emf(
    alkalinity_certified,
    titrant_mass,
    dd.measurement,
    dd.temperature,
    analyte_mass,
    totals,
    k_constants,
    emf0_init=emf0_init,
    pH_min=3,
    pH_max=4,
    titrant_normality=titrant_normality,
    titrant_total_sulfate=titrant_total_sulfate,
)
titrant_molinity = cal["x"][0]
totals = calk.core.add_titrant_totals(
    totals,
    titrant_mass,
    analyte_mass,
    titrant_molinity,
    titrant_total_sulfate=titrant_total_sulfate,
)
sr = calk.core.solve_emf(
    titrant_molinity,
    titrant_mass,
    dd.measurement,
    dd.temperature,
    analyte_mass,
    totals,
    k_constants,
    emf0_init=emf0_init,
    pH_min=3,
    pH_max=4,
    titrant_normality=titrant_normality,
)
print(cal["x"][0])
print(fcal["x"][0])
print(sr.ggr.alkalinity, sr.ggr.emf0)
print(fsr.ggr.alkalinity, fsr.ggr.emf0)
print(sr.emf0, sr.alkalinity)
print(fsr.emf0, fsr.alkalinity)

# %% Visualise alkalinity guess
fig, ax = plt.subplots(dpi=300)
ax.scatter(
    titrant_mass[~ggr.used] * 1e3,
    ggr.gfunc[~ggr.used] * 1e-6,
    c="xkcd:dark",
)
ax.scatter(
    titrant_mass[ggr.used] * 1e3,
    ggr.gfunc[ggr.used] * 1e-6,
    c="xkcd:strawberry",
)
ax.axline(
    (ggr.intercept_x * 1e3, 0),
    slope=ggr.lr.slope * 1e-9,
    c="xkcd:strawberry",
)
ax.axhline(0, c="k", lw=0.8)
ax.axvline(0, c="k", lw=0.8)
ax.grid(alpha=0.2)
ax.set_xlabel("Titrant mass / g")
ax.set_ylabel("Gran estimate / $10^6$")
fig.tight_layout()

fig, ax = plt.subplots(dpi=300)
ax.scatter(
    titrant_mass[~ggr.used] * 1e3,
    ggr.emf0s[~ggr.used],
    c="xkcd:dark",
)
ax.scatter(
    titrant_mass[ggr.used] * 1e3,
    ggr.emf0s[ggr.used],
    c="xkcd:strawberry",
)
ax.axhline(ggr.emf0, c="xkcd:strawberry")
ax.grid(alpha=0.2)
ax.set_xlabel("Titrant mass / g")
ax.set_ylabel(r"EMF$^\mathrm{o}$ estimate / mV")
fig.tight_layout()
