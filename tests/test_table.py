import copy
import numpy as np
import pandas as pd
import calkulate as calk


def test_table_csv():
    """Check that a table can be imported from CSV, calibrated and solved."""
    tt = calk.Dataset("tests/data/titration_table.csv")
    tt.calibrate()
    tt.solve()
    assert ~np.any(np.isnan(tt.table.alkalinity))
    return tt


def test_table_xlsx():
    """Check that a table can be imported from Excel, calibrated and solved."""
    tt = calk.Dataset("tests/data/titration_table.xlsx")
    tt.calibrate()
    tt.solve()
    assert ~np.any(np.isnan(tt.table.alkalinity))


def test_table_customfunc():
    """Check that a table can be imported from CSV, calibrated and solved, with a user-
    specified import function.
    """
    tt = calk.Dataset("tests/data/titration_table.csv", read_func=pd.read_csv)
    tt.calibrate()
    tt.solve()
    assert ~np.any(np.isnan(tt.table.alkalinity))


tt = test_table_csv()
test_table_xlsx()
test_table_customfunc()


def test_bad_filename():
    """Check that a bad filename doesn't break everything."""
    bad_file_table = tt.table.drop(axis=1, labels="alkalinity")
    i = 0
    bad_file_table.loc[i, "file_name"] = "blah.dat"
    bad_file_test = calk.Dataset(bad_file_table)
    bad_file_test.calibrate()
    bad_file_test.solve()
    assert bad_file_test.titrations[i] is None
    assert np.isnan(bad_file_test.table.loc[i].alkalinity)


test_bad_filename()
