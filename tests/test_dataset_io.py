import calkulate as calk, pandas as pd


def test_read_csv():
    """Can a table can be imported from CSV?"""
    ds = calk.read_csv("tests/data/titration_table.csv")
    assert isinstance(ds, pd.DataFrame)
    assert isinstance(ds, calk.Dataset)


def test_read_excel():
    """Can a table can be imported from Excel?"""
    ds = calk.read_excel("tests/data/titration_table.xlsx")
    assert isinstance(ds, pd.DataFrame)
    assert isinstance(ds, calk.Dataset)


def test_read_dbs():
    """Can a table can be imported from a VINDTA .dbs file?"""
    ds = calk.read_dbs(
        "tests/data/vindta_database.dbs", file_path="tests/data/vindta_database/"
    )
    assert isinstance(ds, pd.DataFrame)
    assert isinstance(ds, calk.Dataset)


# ds = calk.read_csv("tests/data/titration_table.csv").calkulate()
# test_read_csv()
# test_read_excel()
# test_read_dbs()
