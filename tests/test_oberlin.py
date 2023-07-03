import calkulate as calk

# Import metadata
ds = calk.read_csv("tests/data/oberlin/metadf.csv")

# Solve for alkalinity
ds.solve()
