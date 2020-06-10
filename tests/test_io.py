import pandas as pd
import numpy as np
import calkulate as calk


titration_table = pd.read_csv("tests/data/titration_table.csv")
t = calk.types.Potentiometric(titration_table.loc[0])


def check_read_dat(i):
    """Import potentiometric titration data from a text file and check they look like
    they should."""
    titration = calk.types.Potentiometric(titration_table.loc[i])
    assert hasattr(titration, "fname")
    assert titration.fname == titration_table.loc[i].fname
    assert hasattr(titration, "analyte")
    assert hasattr(titration, "titrant")
    assert hasattr(titration, "mixture")
    assert hasattr(titration.titrant, "volume")
    assert hasattr(titration.mixture, "emf")
    assert hasattr(titration.mixture, "temperature")
    assert isinstance(titration.fname, str)
    assert isinstance(titration.analyte, calk.types.components.Analyte)
    assert isinstance(titration.titrant, calk.types.components.Titrant)
    assert isinstance(titration.mixture, calk.types.components.Mixture)
    assert isinstance(titration.titrant.volume, np.ndarray)
    assert isinstance(titration.mixture.emf, np.ndarray)
    assert isinstance(titration.mixture.temperature, np.ndarray)
    assert len(np.shape(titration.titrant.volume)) == 1
    assert len(np.shape(titration.mixture.emf)) == 1
    assert len(np.shape(titration.mixture.temperature)) == 1
    assert (
        np.size(titration.titrant.volume)
        == np.size(titration.mixture.emf)
        == np.size(titration.mixture.temperature)
    )
    assert hasattr(titration.analyte, "mass")
    assert hasattr(titration.titrant, "mass")
    assert hasattr(titration.mixture, "dilution_factor")
    assert np.size(titration.mixture.dilution_factor) == np.size(titration.mixture.emf)


def test_read_dat():
    """Apply check_read_dat to every line of the test file."""
    for i in range(len(titration_table.index)):
        check_read_dat(i)


def test_read_write_read_dat():
    """Import potentiometric titration data from a text file, write them to a new text
    file, and re-import.  Check that the values haven't changed.
    """
    titration = calk.types.Potentiometric(titration_table.loc[0])
    titration.write_dat(titration_table.loc[1].fname, mode="w")
    titration_copy = calk.types.Potentiometric(titration_table.loc[1])
    assert np.all(titration.titrant.volume == titration_copy.titrant.volume)
    assert np.all(titration.mixture.emf == titration_copy.mixture.emf)
    assert np.all(titration.mixture.temperature == titration_copy.mixture.temperature)


test_read_dat()
test_read_write_read_dat()
