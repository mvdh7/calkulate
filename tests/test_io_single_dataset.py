import calkulate as calk

tt = calk.read_csv("tests/data/titration_table.csv")
tt = tt.get_dat_files()
dbs = calk.read_dbs(
    "tests/data/vindta_database.dbs", file_path="tests/data/vindta_database/"
).get_dat_files()


def test_table_csv():
    """Can a table can be imported from CSV?"""
    tt = calk.read_csv("tests/data/titration_table.csv")
    tt = tt.get_dat_files()
    assert isinstance(tt, calk.io.Dataset)


def test_table_excel():
    """Can a table can be imported from Excel?"""
    tt = calk.read_excel("tests/data/titration_table.xlsx")
    assert isinstance(tt, calk.io.Dataset)


def test_table_dbs():
    """Can a table can be imported from a VINDTA .dbs file?"""
    dbs = calk.read_dbs(
        "tests/data/vindta_database.dbs", file_path="tests/data/vindta_database/"
    ).get_dat_files()
    assert isinstance(dbs, calk.io.Dataset)


test_table_csv()
test_table_excel()
test_table_dbs()
