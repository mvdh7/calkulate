# %%
import warnings

import numpy as np
import pandas as pd

import calkulate as calk


fname_dbs = "tests/data/vindta_database.dbs"
fpath_dbs = "tests/data/vindta_database/"
dbs = calk.read_dbs(fname_dbs, file_path=fpath_dbs, analyte_volume=97.7)
dbs["alkalinity_certified"] = np.where(dbs.station == 666, 2215, np.nan)
# dbs.calkulate()
# dbs.to_parquet("tests/data/test_dbs_calk_v23_6_2.parquet")


def test_dbs_calkulate():
    """Can we run Calkulate on a dbs file?"""
    assert isinstance(dbs, pd.DataFrame)
    assert isinstance(dbs, calk.Dataset)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        dbs.calkulate(verbose=False)
    assert "alkalinity" in dbs
    assert np.isclose(
        (dbs.alkalinity - dbs.alkalinity_certified).mean(),
        0,
        rtol=0,
        atol=1e-3,
    )


def test_values():
    """Does the current version return the same values as previous versions?"""
    dbs_legacy = calk.calibrate(
        dbs.copy(),
        verbose=False,
        # These kwargs needed for consistency with versions before v23.7
        dilute_totals_for_ks=True,
        double=False,
        gran_logic="legacy",
    )
    dbs_newver = calk.calibrate(dbs.copy(), verbose=False)
    dbs_23_6_2 = pd.read_parquet("tests/data/test_dbs_calk_v23_6_2.parquet")
    dbs_23_7_0 = pd.read_parquet("tests/data/test_dbs_calk_v23_7_0.parquet")
    assert np.allclose(
        dbs_legacy.titrant_molinity,
        dbs_23_6_2.titrant_molinity,
        rtol=0,
        atol=1e-10,
    )
    L = dbs_23_6_2.alkalinity.notnull()
    assert np.allclose(
        dbs_legacy[L].gran_emf0,
        dbs_23_6_2[L].emf0_gran,
        rtol=0,
        atol=1e-7,
    )
    assert np.allclose(
        dbs_legacy[L].gran_alkalinity,
        dbs_23_6_2[L].alkalinity_gran,
        rtol=0,
        atol=1e-5,
    )
    assert np.allclose(
        dbs_legacy[L].emf0,
        dbs_23_6_2[L].emf0,
        rtol=0,
        atol=1e-5,
    )
    assert np.allclose(
        dbs_legacy[L].alkalinity,
        dbs_23_6_2[L].alkalinity,
        rtol=0,
        atol=1e-5,
    )
    assert np.allclose(
        dbs_newver.titrant_molinity,
        dbs_23_7_0.titrant_molinity,
        rtol=0,
        atol=1e-10,
    )
    L = dbs_23_6_2.alkalinity.notnull()
    assert np.allclose(
        dbs_newver[L].gran_emf0,
        dbs_23_7_0[L].gran_emf0,
        rtol=0,
        atol=1e-7,
    )
    assert np.allclose(
        dbs_newver[L].gran_alkalinity,
        dbs_23_7_0[L].gran_alkalinity,
        rtol=0,
        atol=1e-5,
    )
    assert np.allclose(
        dbs_newver[L].emf0,
        dbs_23_7_0[L].emf0,
        rtol=0,
        atol=1e-5,
    )
    assert np.allclose(
        dbs_newver[L].alkalinity,
        dbs_23_7_0[L].alkalinity,
        rtol=0,
        atol=1e-5,
    )


def test_dbs_to_Titration():
    """Can we convert a dbs row to a Titration and get the same results?"""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        dbs.calkulate(verbose=False)
    ix = 20
    tt = dbs.to_Titration(ix)
    assert np.isclose(
        tt.alkalinity, dbs.loc[ix, "alkalinity"], rtol=0, atol=1e-12
    )
    assert np.isclose(tt.emf0, dbs.loc[ix, "emf0"], rtol=0, atol=1e-12)


# test_dbs_calkulate()
# test_dbs_to_Titration()
# test_values()
