import numpy as np
import calkulate as calk


fname = "tests/data/CRM-144-0435-4.dat"  # used by all the tests


def test_read_dat():
    """Import potentiometric titration data from a text file and check they look like
    they should."""
    titration = calk.types.Potentiometric(fname)
    assert hasattr(titration, "titrant")
    assert hasattr(titration, "mixture")
    assert hasattr(titration.titrant, "volume")
    assert hasattr(titration.mixture, "emf")
    assert hasattr(titration.mixture, "temperature")
    assert isinstance(titration.titrant.volume, np.ndarray)
    assert isinstance(titration.mixture.emf, np.ndarray)
    assert isinstance(titration.mixture.temperature, np.ndarray)
    assert len(np.shape(titration.titrant.volume)) == 1
    assert len(np.shape(titration.mixture.emf)) == 1
    assert len(np.shape(titration.mixture.temperature)) == 1
    assert (
        len(titration.titrant.volume)
        == len(titration.mixture.emf)
        == len(titration.mixture.temperature)
    )


def test_read_write_read_dat():
    """Import potentiometric titration data from a text file, write them to a new text
    file, and re-import.  Check that the values haven't changed.
    """
    titration = calk.types.Potentiometric(fname)
    fname_copy = fname.replace(".dat", "-copy.dat")
    titration.write_dat(fname_copy, mode="w")
    titration_copy = calk.types.Potentiometric(fname_copy)
    assert np.all(titration.titrant.volume == titration_copy.titrant.volume)
    assert np.all(titration.mixture.emf == titration_copy.mixture.emf)
    assert np.all(titration.mixture.temperature == titration_copy.mixture.temperature)
