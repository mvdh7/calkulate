import calkulate as calk

ds = calk.read_csv("tests/data/titration_table.csv")


calk.plot.single.emf(ds.loc[0])
