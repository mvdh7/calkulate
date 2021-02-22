import warnings
import calkulate as calk, pandas as pd, numpy as np

fname_dbs = "tests/data/vindta_database.dbs"
fpath_dbs = "tests/data/vindta_database/"
dbs = calk.read_dbs(fname_dbs, file_path=fpath_dbs, analyte_volume=97.7)


def test_dbs_calkulate():
    """Can we run Calkulate on a dbs file?"""
    assert isinstance(dbs, pd.DataFrame)
    assert isinstance(dbs, calk.Dataset)
    dbs["alkalinity_certified"] = np.where(dbs.station == 666, 2215, np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        dbs.calkulate(verbose=False)
    assert "alkalinity" in dbs
    assert np.isclose(
        (dbs.alkalinity - dbs.alkalinity_certified).mean(), 0, rtol=0, atol=1e-3
    )


# test_dbs_calkulate()
