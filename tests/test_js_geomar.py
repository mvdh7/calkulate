# %%
import numpy as np
import pandas as pd

import calkulate as calk


# Import dataset and use dataset functions to calibrate and solve
tf = calk.read_csv("tests/data/js-geomar/js-geomar.csv")
tf.calkulate()
tf = pd.DataFrame(tf)


def test_js_dataset_calkulate():
    assert isinstance(tf, pd.DataFrame)
    assert tf.alkalinity.notnull().all()
    assert tf.emf0.notnull().all()
    crm = tf.alkalinity_certified.notnull()
    assert crm.sum() == 1
    assert np.allclose(tf.alkalinity_certified[crm], tf.alkalinity[crm])


# test_js_dataset_calkulate()
