import numpy as np
import pandas as pd
import calkulate as calk

titration_table = pd.read_csv("tests/data/titration_table.csv")
tt0 = calk.Titration(titration_table.loc[3])
tt1 = calk.Titration(titration_table.loc[4])


def test_Dickson1981():
    # No phosphate
    tt0.solve(pH_range=(0, 10))
    assert np.all(
        np.abs(tt0.solved["alkalinity_points"] * 1e6 - tt0.analyte.alkalinity_certified)
        < 1e-3
    )
    assert tt0.solved["alkalinity_std"] < 1e-9
    # With phosphate - ignoring Dickson's potential typos in a couple of the values
    tt1.solve(pH_range=(0, 10))
    good_Dickson = (np.abs(tt1.titrant.mass - 4.5e-4) > 1e-12) & (
        np.abs(tt1.titrant.mass - 6e-4) > 1e-12
    )
    assert np.all(
        np.abs(
            tt1.solved["alkalinity_points"] * 1e6 - tt1.analyte.alkalinity_certified
        )[good_Dickson]
        < 1e-2
    )
    assert np.std(tt1.solved["alkalinity_points"][good_Dickson]) < 1e-8


test_Dickson1981()
