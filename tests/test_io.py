import pandas as pd
import numpy as np
import calkulate as calk


titration_table = pd.DataFrame(
    {"fname": ["tests/data/CRM-144-0435-4.dat", "tests/data/CRM-144-0435-4-copy.dat"],}
)


def test_read_dat():
    """Import potentiometric titration data from a text file and check they look like
    they should."""
    titration = calk.types.Potentiometric(titration_table.loc[0])
    assert hasattr(titration, "fname")
    assert titration.fname == titration_table.loc[0].fname
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
        len(titration.titrant.volume)
        == len(titration.mixture.emf)
        == len(titration.mixture.temperature)
    )


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
