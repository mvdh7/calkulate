import numpy as np, pandas as pd
from ..solve import complete_emf, complete_pH


def solve(dataset):
    """Determine alkalinity for all samples that have a titrant_molinity value."""
    alkalinity = {}
    emf0 = {}
    for i, row in dataset.iterrows():
        if row.titration is not None and ~np.isnan(row.titrant_molinity):
            if row.measurement_type == "pH":
                solved = complete_pH(row.titration, row)
                alkalinity[i] = solved["x"][0]
                emf0[i] = None
            elif row.measurement_type == "emf":
                solved = complete_emf(row.titration, row)
                alkalinity[i], emf0[i] = solved["x"]
        else:
            alkalinity[i] = None
            emf0[i] = None
    dataset["alkalinity"] = pd.Series(alkalinity) * 1e6
    return dataset
