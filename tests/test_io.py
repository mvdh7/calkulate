import numpy as np
import calkulate as calk


def test_read_dat():
    acid_volume, emf, temperature = calk.io.read_dat("tests/data/CRM-144-0435-4.dat")
    assert isinstance(acid_volume, np.ndarray)
    assert isinstance(emf, np.ndarray)
    assert isinstance(temperature, np.ndarray)
    assert len(np.shape(acid_volume)) == 1
    assert len(np.shape(emf)) == 1
    assert len(np.shape(temperature)) == 1
    assert len(acid_volume) == len(emf) == len(temperature)
