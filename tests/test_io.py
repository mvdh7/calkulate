import pandas as pd
import numpy as np
import calkulate as calk


titration_table = pd.read_csv("tests/data/titration_table.csv")
tt = calk.Titration(titration_table.loc[3])


def check_read_dat(i):
    """Import potentiometric titration data from a text file and check they look like
    they should."""
    titration = calk.Titration(titration_table.loc[i])
    assert isinstance(titration, calk.types.Titration)
    assert hasattr(titration, "fname")
    assert (
        titration.fname
        == titration_table.loc[i].file_path + titration_table.loc[i].file_name
    )
    assert hasattr(titration, "analyte")
    assert hasattr(titration, "titrant")
    assert hasattr(titration, "mixture")
    assert hasattr(titration.mixture, "temperature")
    assert isinstance(titration.fname, str)
    assert isinstance(titration.analyte, calk.types.components.Analyte)
    assert isinstance(titration.titrant, calk.types.components.Titrant)
    assert isinstance(titration.mixture, calk.types.components.Mixture)
    assert isinstance(titration.mixture.temperature, np.ndarray)
    if titration.measurement_type == "EMF":
        assert hasattr(titration.mixture, "emf")
        assert isinstance(titration.mixture.emf, np.ndarray)
        assert len(np.shape(titration.mixture.emf)) == 1
        assert np.size(titration.mixture.emf) == np.size(titration.mixture.temperature)
    assert len(np.shape(titration.mixture.temperature)) == 1
    assert hasattr(titration.analyte, "mass")
    assert hasattr(titration.titrant, "mass")
    assert hasattr(titration.mixture, "dilution_factor")
    assert hasattr(titration.mixture, "total_salts")
    assert hasattr(titration.mixture, "equilibrium_constants")
    assert isinstance(titration.mixture.total_salts, dict)
    assert isinstance(titration.mixture.equilibrium_constants, dict)
    return titration


def check_printing(tt):
    """Check that printing the various classes doesn't throw any errors."""
    print(tt)
    print(tt.analyte)
    print(tt.mixture)
    print(tt.titrant)
    print(tt.settings)


def test_read_dat():
    """Apply check_read_dat to every line of the test file."""
    for i in range(len(titration_table.index)):
        tt = check_read_dat(i)
        check_printing(tt)


def test_read_write_read_dat():
    """Import potentiometric titration data from a text file, write them to a new text
    file, and re-import.  Check that the values haven't changed.
    """
    titration = calk.Titration(titration_table.loc[0])
    titration.write_dat(
        titration_table.loc[1].file_path + titration_table.loc[1].file_name, mode="w"
    )
    titration_copy = calk.Titration(titration_table.loc[1])
    assert np.all(titration.titrant.volume == titration_copy.titrant.volume)
    assert np.all(titration.mixture.emf == titration_copy.mixture.emf)
    assert np.all(titration.mixture.temperature == titration_copy.mixture.temperature)


test_read_dat()
test_read_write_read_dat()
