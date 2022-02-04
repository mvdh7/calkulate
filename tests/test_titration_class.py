import pandas as pd, numpy as np
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
    file_name='0-0  0  (0)calk-3-5.dat',
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
st.get_dic_loss()


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


def test_dic_loss():
    """Does the DIC-loss function work as expected?"""
    assert not hasattr(tt, "k_dic_loss")
    tt.get_dic_loss()
    assert hasattr(tt, "k_dic_loss")
    assert "dic_loss_modelled" in tt.titration
    assert "fCO2_loss_modelled" in tt.titration
    alkalinity0 = tt.alkalinity
    tt.titration["dic"] = tt.titration.dic_loss_modelled
    tt.solve()
    assert not pd.isnull(tt.alkalinity)
    assert not np.isclose(tt.alkalinity, alkalinity0)
    assert hasattr(st, "k_dic_loss")
    assert np.isclose(st.k_dic_loss, 0, rtol=0, atol=1e-10)


def test_plots():
    """Do the plots throw no errors?"""
    tt.plot_emf()
    tt.plot_pH()
    tt.plot_gran_emf0()
    tt.plot_gran_alkalinity()
    tt.plot_alkalinity()
    tt.plot_components()
    tt.plot_components(log_scale=False)
    tt.plot_dic_loss()
    tt.plot_fCO2_loss()


# test_class()
# test_class_simulated()
# test_dic_loss()
# test_plots()
