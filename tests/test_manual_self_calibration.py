import copy
import calkulate as calk, numpy as np

# Import data and self-calibrate
tf = calk.read_csv("tests/data/titration_table.csv").prepare()
tf.calibrate()
tf["titrant_molinity"] = tf["titrant_molinity_here"]  # just for these tests
tf.solve()


def test_manual_self_calibration():
    """Can we successfully calibrate samples with themselves?"""
    l = ~np.isnan(tf.alkalinity)
    assert np.all(
        np.isclose(tf.alkalinity_certified[l], tf.alkalinity[l], rtol=0, atol=1e-6)
    )


test_manual_self_calibration()


def test_Dickson_via_EMF():
    """Can we still return correct alkalinity and EMF from the Dickson data if we first
    convert it from pH to EMF with some arbitrary EMF0?"""
    emf0_i = 678
    tf2 = copy.deepcopy(tf)
    i = 5
    tf2.loc[i, "measurement_type"] = "emf"
    tf2.loc[i].titration["emf"] = calk.convert.h_to_emf(
        10.0 ** -tf2.loc[i].titration["pH"], emf0_i, tf2.loc[i].titration["temperature"]
    )
    alkalinity, emf0 = calk.solvers.complete_emf(tf2.loc[i])["x"]
    assert np.isclose(
        alkalinity * 1e6, tf2.loc[i].alkalinity_certified, rtol=0, atol=1e-3
    )
    assert np.isclose(emf0, emf0_i, rtol=0, atol=1e-4)


test_Dickson_via_EMF()
