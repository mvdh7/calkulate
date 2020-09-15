import calkulate as calk, pandas as pd, numpy as np

file_path = "tests/data/"
file_names = [
    "PC_LIMS_Report-BATCH138-20200317-135120.txt",
    "PC_LIMS_Report-SEA2-20200317-130328.txt",
]


def test_read_TiTouch():
    """Can we import single .dat files into the expected Titration?"""
    dat0 = calk.read_dat(file_path + file_names[0], **calk.kwargs_TiTouch)
    assert isinstance(dat0, pd.DataFrame)
    assert isinstance(dat0, calk.Titration)
    dat1 = calk.read_dat(file_path + file_names[1], **calk.kwargs_TiTouch)
    assert isinstance(dat1, pd.DataFrame)
    assert isinstance(dat1, calk.Titration)


test_read_TiTouch()


def test_TiTouch_dataset():
    """Can we import a Dataset of Ti-Touch-style files and calibrate them?"""
    ds = calk.read_excel("tests/data/TiTouch.xlsx").get_titrations(
        read_dat_kwargs=calk.kwargs_TiTouch
    )
    assert isinstance(ds, calk.Dataset)
    ds = calk.read_excel("tests/data/TiTouch.xlsx").prepare(
        read_dat_kwargs=calk.kwargs_TiTouch
    )
    assert isinstance(ds, calk.Dataset)
    ds = calk.read_excel("tests/data/TiTouch.xlsx").calkulate(
        read_dat_kwargs=calk.kwargs_TiTouch
    )
    assert isinstance(ds, calk.Dataset)
    crm = ~np.isnan(ds.alkalinity_certified)
    print(ds.alkalinity[crm] - ds.alkalinity_certified[crm])
    assert np.all(
        np.isclose(ds.alkalinity[crm], ds.alkalinity_certified[crm], rtol=0, atol=1e-8)
    )
    assert np.all((ds.alkalinity > 1500) & (ds.alkalinity < 3000))
    assert np.all((ds.titrant_molinity > 0.08) & (ds.titrant_molinity < 0.12))
    return ds


ds = test_TiTouch_dataset()
