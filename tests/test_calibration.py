import pandas as pd, numpy as np, calkulate as calk


titration_table = pd.read_csv("tests/data/titration_table.csv")
tt = calk.Titration(titration_table.loc[0])


def test_self_calibration():
    """Test whether solving a sample that you have calibrated with itself returns the
    same certified alkalinity value (within solver tolerance).
    """
    tt.calibrate()
    tt.set_own_titrant_molinity()
    tt.solve()
    assert np.abs(tt.analyte.alkalinity - tt.analyte.alkalinity_certified) < 1e-6


test_self_calibration()
