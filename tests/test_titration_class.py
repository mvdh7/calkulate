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
