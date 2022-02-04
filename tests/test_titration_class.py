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
    # file_name="seawater-CRM-144.dat",
    file_name="0-0  0  (0)calk-3-5.dat",
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

# Run the main function
tt.get_dic_loss()
tt.titration["dic"] = tt.titration.dic_loss_modelled
tt.solve()

#%% Draw figure
k_dic_loss, loss_hires = tt._get_dic_loss_hires()
ttt = tt.titration


from matplotlib import pyplot as plt

fig, axs = plt.subplots(dpi=300, nrows=2, figsize=(6, 5))


ax = axs[0]
ax.plot(
    ttt.titrant_mass * 1000,
    ttt.dic.iloc[0] * 1e6 * ttt.dilution_factor,
    c="k",
    label="Dilution only",
    alpha=0.8,
)
ax.scatter(
    ttt.titrant_mass * 1000,
    ttt.dic_loss,
    s=25,
    label="Calc. from pH",
    c="xkcd:slate",
    alpha=0.8,
    edgecolor="none",
)
ax.fill_between(
    ttt.titrant_mass * 1000,
    ttt.dic_loss_lo,
    y2=ttt.dic_loss_hi,
    data=ttt,
    alpha=0.2,
    label="Calc. uncertainty",
    zorder=-1,
    color="xkcd:slate",
    edgecolor="none",
)
ax.plot(
    loss_hires[loss_hires.pH >= tt.split_pH].titrant_mass * 1000,
    loss_hires[loss_hires.pH >= tt.split_pH].dic,
    label="Model, 'fitted'",
    c="xkcd:teal blue",
    alpha=0.8,
)
ax.plot(
    loss_hires[loss_hires.pH < tt.split_pH].titrant_mass * 1000,
    loss_hires[loss_hires.pH < tt.split_pH].dic,
    label="Model, projected",
    c="xkcd:brownish orange",
    alpha=0.8,
)
ax.set_ylim([loss_hires.dic.min() * 0.96, loss_hires.dic.max() * 1.04])
ax.set_title("$k$(DIC loss) = {:.2f}".format(tt.k_dic_loss))
ax.set_xlabel("Titrant mass / g")
ax.set_ylabel("DIC / $\mu$mol/kg")
ax.legend(fontsize=7, loc="lower left")
ax.grid(alpha=0.3)


ax = axs[1]
ax.scatter(
    ttt.titrant_mass * 1000,
    (ttt.fCO2_loss - tt.fCO2_air) / 1000,
    label="Calc. from pH",
    s=25,
    color="xkcd:slate",
    alpha=0.8,
    edgecolor="none",
)
ax.fill_between(
    ttt.titrant_mass * 1000,
    ttt.fCO2_loss_lo / 1000,
    y2=ttt.fCO2_loss_hi / 1000,
    alpha=0.2,
    label="Calc. uncertainty",
    zorder=-1,
    color="xkcd:slate",
    edgecolor="none",
)
ax.plot(
    loss_hires[loss_hires.pH >= tt.split_pH].titrant_mass * 1000,
    loss_hires[loss_hires.pH >= tt.split_pH].delta_fCO2 / 1000,
    label="Model, 'fitted'",
    c="xkcd:teal blue",
    alpha=0.8,
)
ax.plot(
    loss_hires[loss_hires.pH < tt.split_pH].titrant_mass * 1000,
    loss_hires[loss_hires.pH < tt.split_pH].delta_fCO2 / 1000,
    label="Model, projected",
    c="xkcd:brownish orange",
    alpha=0.8,
)
ax.set_ylim([0, loss_hires.delta_fCO2.max() * 1.5e-3])
ax.set_xlabel("Titrant mass / g")
ax.set_ylabel("$\Delta f$CO$_2$ / matm")
ax.legend(fontsize=7, loc="upper left")
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("tests/figures/test_dic_loss.png")

#%%
# fig, ax = plt.subplots(dpi=300)
# ttt.plot.scatter('pH', 'dic_loss', ax=ax)
# ax.set_ylim([1800, 2200])
