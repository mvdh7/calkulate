# %%
import calkulate as calk


def test_irlon():
    tf = calk.read_csv("tests/data/irlon/TA_metadata_20210623.csv")
    tf.calkulate(read_dat_kwargs={"encoding": "unicode_escape"}, verbose=True)
    assert tf.alkalinity.notnull().all()


# test_irlon()
