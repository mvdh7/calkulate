import calkulate as calk, numpy as np


dat = calk.read_dat("tests/data/seawater-CRM-144.dat")
dat_vars = ["titrant_amount", "mixture_measurement", "mixture_temperature"]


def test_read_dat():
    """Can we import a single .dat file into the expected dict?"""
    assert isinstance(dat, dict)
    for var in dat_vars:
        assert var in dat.keys()
        assert isinstance(dat[var], np.ndarray)
        assert dat[var].ndim == 1
        assert dat[var].shape == dat["titrant_amount"].shape


def test_write_reimport_dat_function():
    """Can we write a .dat file from a dict and re-import it and it stays the same,
    using the write_dat function?
    """
    fname_copy = "tests/data/seawater-CRM-144-function.dat"
    calk.write_dat(dat, fname_copy, mode="w")
    dat_copy = calk.read_dat(fname_copy)
    for var in dat_vars:
        assert np.all(np.isclose(dat[var], dat_copy[var]))


def test_write_reimport_dat_method():
    """Can we write a .dat file from a dict and re-import it and it stays the same,
    using the dat.write method?
    """
    fname_copy = "tests/data/seawater-CRM-144-method.dat"
    dat.write(fname_copy, mode="w")
    dat_copy = calk.read_dat(fname_copy)
    for var in dat_vars:
        assert np.all(np.isclose(dat[var], dat_copy[var]))


test_read_dat()
test_write_reimport_dat_function()
test_write_reimport_dat_method()
