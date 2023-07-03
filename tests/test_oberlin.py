import calkulate as calk

# Import metadata
ds = calk.read_csv("tests/data/oberlin/metadf.csv")

# Solve for alkalinity
ds.solve()


def test_oberlin():
    assert "alkalinity" in ds
    assert not ds.alkalinity.isnull().any()


# test_oberlin()
