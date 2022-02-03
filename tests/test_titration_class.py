import pandas as pd
import calkulate as calk

tt_alkalinity = 2238.6
tt_kwargs = dict(
    dic=2031.53,
    total_silicate=2.5,
    total_phosphate=0.31,
)
tt = calk.Titration(
    analyte_mass=0.1,
    salinity=33.571,
    file_name="seawater-CRM-144.dat",
    # file_name="0-0  0  (0)calk-3-5.dat",
    file_path="tests/data/",
    file_prepare_kwargs=tt_kwargs,
)
tt.calkulate(tt_alkalinity)

# Simulate a similar titration
st = calk.Titration(
    analyte_mass=0.1,
    salinity=33.571,
    simulate_alkalinity=tt_alkalinity,
    simulate_kwargs=dict(
        emf0=tt.emf0,
        titrant_molinity=tt.titrant_molinity,
        titrant_mass_stop=4.2e-3,
        titrant_mass_step=0.15e-3,
        **tt_kwargs,
    ),
)
st.calkulate(tt_alkalinity)

# Simulate again using the more user-friendly function
st2 = calk.simulate.titration(
    tt_alkalinity,
    analyte_mass=0.1,
    emf0=tt.emf0,
    titrant_molinity=tt.titrant_molinity,
    titrant_mass_stop=4.2e-3,
    titrant_mass_step=0.15e-3,
    **tt_kwargs,
)


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


def test_class_simulated():
    """Does the simulated Titration object have the expected attributes?"""
    assert hasattr(st, "alkalinity")
    assert hasattr(st, "emf0")
    assert hasattr(st, "alkalinity_gran")
    assert hasattr(st, "emf0_gran")
    assert hasattr(st, "titration")
    assert isinstance(st.alkalinity, float)
    assert isinstance(st.emf0, float)
    assert isinstance(st.alkalinity_gran, float)
    assert isinstance(st.emf0_gran, float)
    assert isinstance(st.titration, pd.DataFrame)


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
# test_class_simulated()
# test_plots()

#%%
from matplotlib import pyplot as plt
import numpy as np
from scipy import interpolate, optimize


def dic_loss_model_fitted(
    k_dic_loss, i_titrant_mass, i_delta_fCO2, dic_start, analyte_mass, step_tm
):
    """Calculate DIC loss for the fitted section of the titration."""
    i_dic = np.full_like(i_delta_fCO2, np.nan)
    i_dic[0] = dic_start * 1e6
    for i in range(len(i_dic) - 1):
        dilute = calk.convert.get_dilution_factor(
            step_tm, analyte_mass + i_titrant_mass[i]
        )
        i_dic[i + 1] = i_dic[i] * dilute - k_dic_loss * i_delta_fCO2[i] * step_tm
    return i_dic


def _lsqfun_dic_loss_model(
    k_dic_loss,
    i_titrant_mass,
    i_delta_fCO2,
    dic_start,
    analyte_mass,
    step_tm,
    i_dic_loss,
):
    return (
        dic_loss_model_fitted(
            k_dic_loss, i_titrant_mass, i_delta_fCO2, dic_start, analyte_mass, step_tm
        )
        - i_dic_loss
    )


def get_fCO2_from_dic_pH(dic, pH, k0, k1, k2):
    """Calculate fCO2 from DIC and pH."""
    h = 10 ** -pH
    CO2aq = dic / (1 + k1 / h + k1 * k2 / h ** 2)
    fCO2 = CO2aq / k0
    return fCO2


def dic_loss_model_future(
    k_dic_loss,
    f_titrant_mass,
    f_pH,
    f_k0,
    f_k1,
    f_k2,
    delta_fCO2_start,
    dic_start,
    pH_start,
    k0_start,
    k1_start,
    k2_start,
    analyte_mass,
    step_tm,
):
    """Forecast future DIC loss, starting from the end of the fitted section."""
    # Prepare empty arrays for model
    f_delta_fCO2 = np.full_like(f_titrant_mass, np.nan)
    f_dic = np.full_like(f_titrant_mass, np.nan)
    # Get first values in the forecast arrays
    dilute = calk.convert.get_dilution_factor(
        step_tm, analyte_mass + f_titrant_mass[0] - step_tm
    )
    f_dic[0] = dic_start * dilute - k_dic_loss * delta_fCO2_start * step_tm
    f_delta_fCO2[0] = (
        get_fCO2_from_dic_pH(f_dic[0], pH_start, k0_start, k1_start, k2_start)
        - fCO2_air
    )
    # Forecast future DIC loss
    for i in range(len(f_dic) - 1):
        dilute = calk.convert.get_dilution_factor(
            step_tm, analyte_mass + f_titrant_mass[i]
        )
        f_dic[i + 1] = f_dic[i] * dilute - k_dic_loss * f_delta_fCO2[i] * step_tm
        f_delta_fCO2[i + 1] = (
            get_fCO2_from_dic_pH(
                f_dic[i + 1], f_pH[i + 1], f_k0[i + 1], f_k1[i + 1], f_k2[i + 1]
            )
            - fCO2_air
        )
    return f_dic, f_delta_fCO2


def get_dic_loss_hires(tt, fCO2_air=450, split_pH=5.5):
    """Fit and forecast high-resolution DIC loss model."""
    # Relabel for convenience
    ttt = tt.titration
    analyte_mass = tt.analyte_mass
    titrant_mass = ttt.titrant_mass.to_numpy()
    dic_loss = ttt.dic_loss.to_numpy()
    pH = ttt.pH.to_numpy()
    dic_start = ttt.dic.iloc[0]
    fCO2_loss = ttt.fCO2_loss.to_numpy()
    k_CO2 = ttt.k_CO2.to_numpy()
    k_carbonic_1 = ttt.k_carbonic_1.to_numpy()
    k_carbonic_2 = ttt.k_carbonic_2.to_numpy()
    # Get delta-fCO2
    delta_fCO2_loss = fCO2_loss - fCO2_air
    # Use titrant_mass as proxy for titration time: generate high-resolution arrays
    step_tm = 1e-6
    a_titrant_mass = np.arange(0, np.max(titrant_mass), step_tm)
    a_pH = interpolate.pchip_interpolate(titrant_mass, pH, a_titrant_mass)
    i_titrant_mass = a_titrant_mass[a_pH >= split_pH]
    f_titrant_mass = a_titrant_mass[a_pH < split_pH]
    # Interpolate other properties to the high-resolution titrant_mass
    i_dic_loss = interpolate.pchip_interpolate(titrant_mass, dic_loss, i_titrant_mass)
    i_delta_fCO2 = interpolate.pchip_interpolate(
        titrant_mass, delta_fCO2_loss, i_titrant_mass
    )
    f_pH = interpolate.pchip_interpolate(titrant_mass, pH, f_titrant_mass)
    f_k0 = interpolate.pchip_interpolate(titrant_mass, k_CO2, f_titrant_mass)
    f_k1 = interpolate.pchip_interpolate(titrant_mass, k_carbonic_1, f_titrant_mass)
    f_k2 = interpolate.pchip_interpolate(titrant_mass, k_carbonic_2, f_titrant_mass)
    # Get mid-way start-points for forecasting
    pH_start = interpolate.pchip_interpolate(titrant_mass, pH, i_titrant_mass[-1])
    k0_start = interpolate.pchip_interpolate(titrant_mass, k_CO2, i_titrant_mass[-1])
    k1_start = interpolate.pchip_interpolate(
        titrant_mass, k_carbonic_1, i_titrant_mass[-1]
    )
    k2_start = interpolate.pchip_interpolate(
        titrant_mass, k_carbonic_2, i_titrant_mass[-1]
    )
    # Find best-fit k_dic_loss
    k_dic_loss_opt_result = optimize.least_squares(
        _lsqfun_dic_loss_model,
        1.0,
        args=(
            i_titrant_mass,
            i_delta_fCO2,
            dic_start,
            analyte_mass,
            step_tm,
            i_dic_loss,
        ),
    )
    k_dic_loss = k_dic_loss_opt_result["x"]
    # Calculate DIC from k_dic_loss in the fitted region
    i_dic = dic_loss_model_fitted(
        k_dic_loss, i_titrant_mass, i_delta_fCO2, dic_start, analyte_mass, step_tm
    )
    # Forecast future DIC loss
    f_dic, f_delta_fCO2 = dic_loss_model_future(
        k_dic_loss,
        f_titrant_mass,
        f_pH,
        f_k0,
        f_k1,
        f_k2,
        i_delta_fCO2[-1],
        i_dic[-1],
        pH_start,
        k0_start,
        k1_start,
        k2_start,
        analyte_mass,
        step_tm,
    )
    return k_dic_loss[0], {
        "titrant_mass": a_titrant_mass,
        "pH": a_pH,
        "dic": np.concatenate((i_dic, f_dic)),
        "delta_fCO2": np.concatenate((i_delta_fCO2, f_delta_fCO2)),
    }


def get_dic_loss(tt, fCO2_air=450, split_pH=5.5):
    """Get final DIC loss values at the titration points to go in the titration df."""
    tt.k_dic_loss, loss_hires = get_dic_loss_hires(
        tt, fCO2_air=fCO2_air, split_pH=split_pH
    )
    tt.titration["dic_loss_modelled"] = interpolate.pchip_interpolate(
        loss_hires["titrant_mass"],
        loss_hires["dic"],
        tt.titration.titrant_mass.to_numpy(),
    )
    tt.titration["fCO2_loss_modelled"] = (
        interpolate.pchip_interpolate(
            loss_hires["titrant_mass"],
            loss_hires["delta_fCO2"],
            tt.titration.titrant_mass.to_numpy(),
        )
        + fCO2_air
    )
    tt.titration["loss_fitted"] = tt.titration.pH >= split_pH


# Main function inputs
fCO2_air = 450
split_pH = 5.5  # fit DIC-loss model only above this pH, forecast at lower pH
k_dic_loss, loss_hires = get_dic_loss_hires(tt, fCO2_air=fCO2_air, split_pH=split_pH)
get_dic_loss(tt, fCO2_air=fCO2_air, split_pH=split_pH)

#%% Draw figure
loss_hires = pd.DataFrame(loss_hires)
ttt = tt.titration
fig, axs = plt.subplots(dpi=300, nrows=2, figsize=(6, 5))


ax = axs[0]
ax.plot(ttt.titrant_mass, ttt.dic * 1e6, c="k", label="Dilution only")
ax.scatter(
    "titrant_mass", "dic_loss", data=ttt, s=20, label="Calc. from pH", c="xkcd:slate"
)
ax.fill_between(
    "titrant_mass",
    "dic_loss_lo",
    y2="dic_loss_hi",
    data=ttt,
    alpha=0.3,
    label="Calc. uncertainty",
    zorder=-1,
    color="xkcd:slate",
    edgecolor="none",
)
# ax.scatter(0, tt.dic)
ax.plot(
    "titrant_mass",
    "dic",
    data=loss_hires[loss_hires.pH >= split_pH],
    label="Model, 'fitted'",
    c="xkcd:teal blue",
)
ax.plot(
    "titrant_mass",
    "dic",
    data=loss_hires[loss_hires.pH < split_pH],
    label="Model, projected",
    c="xkcd:brownish orange",
)
ax.set_ylabel("DIC / $\mu$mol/kg")
ax.legend(fontsize=7)
ax.set_ylim([1450, 2050])
ax.set_title("$k$(DIC loss) = {:.2f}".format(tt.k_dic_loss))

ax = axs[1]
ax.scatter(
    ttt.titrant_mass,
    ttt.fCO2_loss - fCO2_air,
    label="Calc. from pH",
    s=20,
    color="xkcd:slate",
)
ax.fill_between(
    "titrant_mass",
    "fCO2_loss_lo",
    y2="fCO2_loss_hi",
    data=ttt,
    alpha=0.3,
    label="Calc. uncertainty",
    zorder=-1,
    color="xkcd:slate",
    edgecolor="none",
)
ax.plot(
    "titrant_mass",
    "delta_fCO2",
    data=loss_hires[loss_hires.pH >= split_pH],
    label="Model, 'fitted'",
    c="xkcd:teal blue",
)
ax.plot(
    "titrant_mass",
    "delta_fCO2",
    data=loss_hires[loss_hires.pH < split_pH],
    label="Model, projected",
    c="xkcd:brownish orange",
)
ax.set_ylabel("$\Delta f$CO$_2$ / $\mu$atm")
ax.set_ylim([0, 100000])
ax.legend(fontsize=7)
for ax in axs:
    ax.grid(alpha=0.3)
    ax.set_xlabel("Titrant mass / kg")
plt.tight_layout()
plt.savefig("tests/figures/test_dic_loss.png")

#%%
# fig, ax = plt.subplots(dpi=300)
# ttt.plot.scatter('pH', 'dic_loss', ax=ax)
# ax.set_ylim([1800, 2200])
