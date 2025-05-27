import calkulate as calk

ds = {}
for i in [1, 2, 3]:
    ds[i] = calk.read_csv("tests/data/titration_table.csv")
    ds[i]["opt_pH_scale"] = i
    ds[i].calkulate()

# solving on different pH scales doesn't work yet
