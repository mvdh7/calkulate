import pandas as pd
import calkulate as calk

tt = calk.Titration(
    file_name="seawater-CRM-144.dat",
    file_path="tests/data/",
    analyte_mass=0.1,
    dic=2031.53,
    total_silicate=2.5,
    total_phosphate=0.31,
)
tt.calibrate(2238.6)


def test_class():
    """Does the Titration object have the expected attributes?"""
    assert hasattr(tt, "alkalinity")
    assert hasattr(tt, "emf0")
    assert hasattr(tt, "alkalinity_gran")
    assert hasattr(tt, "emf0_gran")
    assert hasattr(tt, "titration")
    assert isinstance(tt.alkalinity, float)
    assert isinstance(tt.emf0, float)
    assert isinstance(tt.alkalinity_gran, float)
    assert isinstance(tt.emf0_gran, float)
    assert isinstance(tt.titration, pd.DataFrame)


def test_plots():
    """Do the plots throw no errors?"""
    tt.plot_emf()
    tt.plot_pH()
    tt.plot_gran_emf0()
    tt.plot_gran_alkalinity()
    tt.plot_alkalinity()
    tt.plot_components()
    tt.plot_components(log_scale=False)


# test_class()
# test_plots()

#%%
from matplotlib import pyplot as plt
import numpy as np
from scipy import interpolate

ttt = tt.titration
fCO2_air = 450
ttt["delta_fCO2_loss"] = ttt.fCO2_loss - fCO2_air

# Titrant mass as proxy for titration time
dx = 1e-6
ix = np.arange(0, 0.0018, dx)
idilute = calk.convert.get_dilution_factor(dx, tt.analyte_mass)
idic_raw = interpolate.pchip_interpolate(
    ttt.titrant_mass.values, ttt.dic_loss.values, ix
)
# Interpolated delta-fCO2
iy = interpolate.pchip_interpolate(
    ttt.titrant_mass.values, ttt.delta_fCO2_loss.values, ix
)

# Model loss of DIC
k_loss = 2.9e-6
idic = np.full_like(iy, np.nan)
idic[0] = ttt.dic_loss.iloc[0]
idic[0] = ttt.dic.iloc[0] * 1e6
for i in range(1, len(idic)):
    idic[i] = idic[i - 1] * idilute - k_loss * iy[i - 1]

# Forecast future DIC loss
def get_delta_fCO2_from_dic_pH(dic, pH, k0, k1, k2):
    h = 10 ** -pH
    CO2aq = dic / (1 + k1 / h + k1 * k2 / h ** 2)
    fCO2 = CO2aq / k0
    return fCO2 - fCO2_air


fx = np.arange(0.0018, 0.0043, dx)
fdic = np.full_like(fx, np.nan)
fdic[0] = idic[-1]

fpH = interpolate.pchip_interpolate(ttt.titrant_mass.values, ttt.pH.values, fx)
fk0 = interpolate.pchip_interpolate(ttt.titrant_mass.values, ttt.k_CO2.values, fx)
fk1 = interpolate.pchip_interpolate(
    ttt.titrant_mass.values, ttt.k_carbonic_1.values, fx
)
fk2 = interpolate.pchip_interpolate(
    ttt.titrant_mass.values, ttt.k_carbonic_2.values, fx
)

fy = get_delta_fCO2_from_dic_pH(fdic, fpH, fk0, fk1, fk2)


for i in range(len(fdic) - 1):
    fdic[i + 1] = fdic[i] * idilute - k_loss * fy[i]
    fy[i + 1] = get_delta_fCO2_from_dic_pH(
        fdic[i + 1], fpH[i + 1], fk0[i + 1], fk1[i + 1], fk2[i + 1]
    )

# Draw figure
fig, axs = plt.subplots(dpi=300, nrows=2)
ax = axs[0]
ax.fill_between("titrant_mass", "dic_loss_lo", y2="dic_loss_hi", data=ttt, alpha=0.3)
ax.plot(ttt.titrant_mass, ttt.dic * 1e6, c="k")
ax.plot("titrant_mass", "dic_loss", data=ttt)
ax.scatter(0, tt.titration.dic.iloc[0] * 1e6)
ax.plot(ix, idic)
ax.plot(fx, fdic)
ax.set_ylim([800, 2500])
ax.set_ylim([1800, 2100])
# ax.set_ylim([1600, 2100])
ax = axs[1]
ax.plot("titrant_mass", "delta_fCO2_loss", data=ttt)
ax.plot(ix, iy)
ax.plot(fx, fy)
ax.set_ylim([0, 100000])
for ax in axs:
    ax.grid(alpha=0.3)
plt.tight_layout()

#%%
