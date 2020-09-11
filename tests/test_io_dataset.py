import calkulate as calk, pandas as pd


def test_read_csv():
    """Can a table can be imported from CSV?"""
    tf = calk.read_csv("tests/data/titration_table.csv")
    assert isinstance(tf, pd.DataFrame)
    assert isinstance(tf, calk.Dataset)


def test_read_excel():
    """Can a table can be imported from Excel?"""
    tf = calk.read_excel("tests/data/titration_table.xlsx")
    assert isinstance(tf, pd.DataFrame)
    assert isinstance(tf, calk.Dataset)


def test_read_dbs():
    """Can a table can be imported from a VINDTA .dbs file?"""
    tf = calk.read_dbs(
        "tests/data/vindta_database.dbs", file_path="tests/data/vindta_database/"
    )
    assert isinstance(tf, pd.DataFrame)
    assert isinstance(tf, calk.Dataset)


def test_get_titrations_method():
    """Can we import titrations with the method?"""
    tf = calk.read_csv("tests/data/titration_table.csv")
    tf.get_titrations()
    assert "titration" in tf
    assert isinstance(tf.titration[0], calk.Titration)


def test_get_titrations_function():
    """Can we import titrations with the function?"""
    tf = calk.read_csv("tests/data/titration_table.csv")
    tf = calk.get_titrations(tf)
    assert "titration" in tf
    assert isinstance(tf.titration[0], calk.Titration)


test_read_csv()
test_read_excel()
test_read_dbs()
test_get_titrations_method()
test_get_titrations_function()
