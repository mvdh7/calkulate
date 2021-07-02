import warnings
import calkulate as calk, pandas as pd, numpy as np

fname_dbs = "tests/data/vindta_database.dbs"
fpath_dbs = "tests/data/vindta_database/"
dbs = calk.read_dbs(fname_dbs, file_path=fpath_dbs, analyte_volume=97.7)
dbs["alkalinity_certified"] = np.where(dbs.station == 666, 2215, np.nan)


def test_dbs_calkulate():
    """Can we run Calkulate on a dbs file?"""
    assert isinstance(dbs, pd.DataFrame)
    assert isinstance(dbs, calk.Dataset)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        dbs.calkulate(verbose=False)
    assert "alkalinity" in dbs
    assert np.isclose(
        (dbs.alkalinity - dbs.alkalinity_certified).mean(), 0, rtol=0, atol=1e-3
    )


def test_dbs_to_Titration():
    """Can we convert a dbs row to a Titration and get the same results?"""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        dbs.calkulate(verbose=False)
    ix = 20
    tt = dbs.to_Titration(ix)
    assert tt.solved
    assert not tt.calibrated
    assert np.isclose(tt.alkalinity, dbs.loc[ix, "alkalinity"], rtol=0, atol=1e-12)
    assert np.isclose(tt.emf0, dbs.loc[ix, "emf0"], rtol=0, atol=1e-12)
    ix_crm = 88  # must be a CRM!
    tt_crm = dbs.to_Titration(ix_crm)
    assert not tt_crm.calibrated
    tt_crm.calibrate(tt_crm.alkalinity_certified)
    assert tt_crm.calibrated
    assert tt_crm.solved
    assert np.isclose(
        tt_crm.titrant_molinity,
        dbs.loc[ix_crm, "titrant_molinity_here"],
        rtol=0,
        atol=1e-12,
    )


# test_dbs_calkulate()
# test_dbs_to_Titration()
