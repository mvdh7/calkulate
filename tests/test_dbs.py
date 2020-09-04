import calkulate as calk, pandas as pd, numpy as np

fname_dbs = "tests/data/vindta_database.dbs"
fpath_dbs = "tests/data/vindta_database/"


def test_dbs_import():
    # Import the dbs
    dbs = calk.io.read_dbs(fname_dbs, file_path=fpath_dbs, analyte_volume=97.7)
    assert isinstance(dbs, pd.DataFrame)
    # Add extra columns
    dbs["alkalinity_certified"] = np.where(dbs.station == 666, 2215, np.nan)
    # Convert to titration table
    tdata = calk.Dataset(dbs)
    assert isinstance(tdata, calk.Dataset)
    assert isinstance(tdata.table, pd.DataFrame)
    assert isinstance(tdata.titrations, dict)
    tdata.calibrate_and_solve()
    tdata.plot_calibration()
    return tdata


def test_dbs_direct():
    tdata = calk.Dataset(fname_dbs, file_path=fpath_dbs)
    assert isinstance(tdata, calk.Dataset)
    assert isinstance(tdata.table, pd.DataFrame)
    assert isinstance(tdata.titrations, dict)
    return tdata


test_dbs_direct()
tdata = test_dbs_import()
