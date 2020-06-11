import pandas as pd
import numpy as np
import calkulate as calk


titration_table = pd.read_csv("tests/data/titration_table.csv")
tt = calk.Potentiometric(titration_table.loc[0])


def test_self_calibration():
    """Test whether solving a sample that you have calibrated with itself returns the
    same certified alkalinity value (within solver tolerance).
    """
    tt.calibrate()
    tt.set_own_titrant_molinity()
    tt.solve()
    assert np.abs(tt.analyte.alkalinity - tt.analyte.alkalinity_certified) < 1e-7


test_self_calibration()
