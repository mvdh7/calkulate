import calkulate as calk, pandas as pd, numpy as np


filestem = "tests/data/seawater-CRM-144"
dat = calk.read_dat(filestem + ".dat")
dat_vars = ["titrant_amount", "measurement", "temperature"]


def test_read_dat():
    """Can we import a single .dat file into the expected dict?"""
    assert isinstance(dat, pd.DataFrame)
    assert isinstance(dat, calk.Titration)
    for var in dat_vars:
        assert var in dat.columns
        assert dat[var].dtype == float


def test_write_reimport_dat_function():
    """Can we write a .dat file from a dict and re-import it and it stays the same,
    using the write_dat function?
    """
    fname_copy = filestem + "-function.dat"
    calk.to_dat(dat, fname_copy, mode="w")
    dat_copy = calk.read_dat(fname_copy)
    for var in dat_vars:
        assert np.all(np.isclose(dat[var], dat_copy[var]))


def test_write_reimport_dat_method():
    """Can we write a .dat file from a dict and re-import it and it stays the same,
    using the dat.write method?
    """
    fname_copy = filestem + "-method.dat"
    dat.to_dat(fname_copy, mode="w")
    dat_copy = calk.read_dat(fname_copy)
    for var in dat_vars:
        assert np.all(np.isclose(dat[var], dat_copy[var]))


test_read_dat()
test_write_reimport_dat_function()
test_write_reimport_dat_method()
