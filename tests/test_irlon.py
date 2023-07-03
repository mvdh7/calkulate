import calkulate as calk


tf = calk.read_csv("tests/data/irlon/TA_metadata_20210623.csv")
tf.calkulate(read_dat_kwargs={"encoding": "unicode_escape"}, verbose=True)


def test_irlon():
    assert not tf.alkalinity.isnull().any()


# test_irlon()
