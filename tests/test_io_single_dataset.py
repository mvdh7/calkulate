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
    assert isinstance(tt, calk.datasets.Dataset)
    assert "dat_dict" in tt


def test_get_dat_twice():
    """Does re-doing get_dat_files() throw no errors?"""
    tt = calk.read_csv("tests/data/titration_table.csv")
    tt = tt.get_dat_files()
    tt = tt.get_dat_files()
    tt = tt.get_dat_files()


def test_get_dat_inplace():
    """Does the in-place get_dat_files() version work?"""
    tt = calk.read_csv("tests/data/titration_table.csv").get_dat_files()
    tt.drop("dat_dict", axis=1, inplace=True)
    assert "dat_dict" not in tt
    tt.get_dat_files()
    assert "dat_dict" in tt


def test_table_excel():
    """Can a table can be imported from Excel?"""
    tt = calk.read_excel("tests/data/titration_table.xlsx")
    assert isinstance(tt, calk.datasets.Dataset)


def test_table_dbs():
    """Can a table can be imported from a VINDTA .dbs file?"""
    dbs = calk.read_dbs(
        "tests/data/vindta_database.dbs", file_path="tests/data/vindta_database/"
    ).get_dat_files()
    assert isinstance(dbs, calk.datasets.Dataset)


test_table_csv()
test_get_dat_twice()
test_get_dat_inplace()
test_table_excel()
test_table_dbs()
