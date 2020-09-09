import calkulate as calk

tt = calk.read_csv("tests/data/titration_table.csv")
tt = tt.get_dat_files()


def test_table_csv():
    """Can a table can be imported from CSV, calibrated and solved?"""
    tt = calk.read_csv("tests/data/titration_table.csv")
    tt = tt.get_dat_files()
    assert isinstance(tt, calk.io.Dataset)


def test_table_excel():
    """Check that a table can be imported from Excel, calibrated and solved."""
    tt = calk.read_excel("tests/data/titration_table.xlsx")
    assert isinstance(tt, calk.io.Dataset)


# def test_table_customfunc():
#     """Check that a table can be imported from CSV, calibrated and solved, with a user-
#     specified import function.
#     """
#     tt = calk.io.Dataset("tests/data/titration_table.csv", read_func=pd.read_csv)
#     # tt.calibrate()
#     # tt.solve()
#     # assert ~np.any(np.isnan(tt.table.alkalinity))


test_table_csv()
test_table_excel()
