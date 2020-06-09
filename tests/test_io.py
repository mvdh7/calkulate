import numpy as np
import calkulate as calk


fname = "tests/data/CRM-144-0435-4.dat"


def test_read_dat():
    """Import data from a text file and check the outputs look like they should."""
    acid_volume, emf, temperature = calk.io.read_dat(fname)
    assert isinstance(acid_volume, np.ndarray)
    assert isinstance(emf, np.ndarray)
    assert isinstance(temperature, np.ndarray)
    assert len(np.shape(acid_volume)) == 1
    assert len(np.shape(emf)) == 1
    assert len(np.shape(temperature)) == 1
    assert len(acid_volume) == len(emf) == len(temperature)


def test_read_write_read_dat():
    """Import data from a text file, write it back to a text file, and re-import.
    Check that the values haven't changed.
    """
    read_vars = calk.io.read_dat(fname)
    fname_copy = fname.replace(".dat", "-copy.dat")
    calk.io.write_dat(fname_copy, *read_vars, mode="w")
    read_vars_again = calk.io.read_dat(fname_copy)
    for i in range(len(read_vars)):
        assert np.all(read_vars[i] == read_vars_again[i])
