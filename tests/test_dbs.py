import calkulate as calk, pandas as pd


def test_dbs_import():
    dbs = calk.io.read_dbs("tests/data/vindta_database.dbs")
    assert isinstance(dbs, pd.DataFrame)
    return dbs


dbs = test_dbs_import()
