# %%
import calkulate as calk


def test_oberlin():
    # Import metadata & solve for alkalinity
    ds = calk.read_csv("tests/data/oberlin/metadf.csv")
    ds.solve()
    assert "alkalinity" in ds
    assert not ds.alkalinity.isnull().any()


# test_oberlin()
