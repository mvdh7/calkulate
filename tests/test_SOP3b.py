import calkulate as calk, numpy as np

file_path = "tests/data/SOP3b/"
ds = calk.read_csv(file_path + "metadata_DicksonSOP3b.csv")
ds["file_path"] = file_path
ds["total_borate"] = 0
ds["opt_k_bisulfate"] = 1
ds["opt_k_fluoride"] = 1
ds.solve(read_dat_kwargs={"delimiter": ","}, pH_range=(0, 14))


def test_SOP3b():
    """Does Calkulate return values close to those reported in the SOP?

    Not exactly, but:
      * the alkalinity is very close (great), but the EMF0 is a bit worse (who cares?);
      * the value of the alkalinity can be shifted around by about ±0.5 μmol/kg by
        switching to different PyCO2SYS options for K-bisulfate and K-fluoride, and it's
        not clear which of these were used in the SOP calculation.
    Given that Calkulate *does* perfectly reproduce the simulated titration data from
    Dickson (1981), where exact values of all equilibrium constants etc. are known, we
    therefore consider the agreement here to be good enough.
    """
    assert np.isclose(ds.alkalinity, 2260.06, rtol=0, atol=0.1)
    assert np.isclose(ds.emf0, 394.401, rtol=0, atol=10)


# test_SOP3b()
