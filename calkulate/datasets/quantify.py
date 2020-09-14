import numpy as np, pandas as pd
from .. import options, solve


def calibrate_all(dataset, pH_range=(3, 4), verbose=options.verbose):
    """Determine titrant molinity for all samples that have a certified alkalinity."""
    titrant_molinity = {}
    for i, row in dataset.iterrows():
        if row.titration is not None and ~np.isnan(row.alkalinity_certified):
            if verbose:
                print("Calkulate: calibrating {}...".format(row.file_name))
            if row.measurement_type == "pH":
                titrant_molinity[i] = solve.calibrate(
                    row, pH_range=pH_range, solver=solve.complete_pH
                )["x"][0]
            else:  # if row.measurement_type == "emf":
                titrant_molinity[i] = solve.calibrate(
                    row, pH_range=pH_range, solver=solve.complete_emf
                )["x"][0]
        else:
            titrant_molinity[i] = None
    dataset["titrant_molinity_here"] = pd.Series(titrant_molinity)
    if verbose:
        print("Calkulate: calibrations complete!")
    return dataset


def solve_all(dataset, pH_range=(3, 4), verbose=options.verbose):
    """Determine alkalinity for all samples that have a titrant molinity value."""
    alkalinity = {}
    emf0 = {}
    for i, row in dataset.iterrows():
        if row.titration is not None and ~np.isnan(row.titrant_molinity):
            if verbose:
                print("Calkulate: solving {}...".format(row.file_name))
            if row.measurement_type == "pH":
                solved = solve.complete_pH(row, pH_range=pH_range)
                alkalinity[i] = solved["x"][0]
                emf0[i] = None
            elif row.measurement_type == "emf":
                solved = solve.complete_emf(row, pH_range=pH_range)
                alkalinity[i], emf0[i] = solved["x"]
        else:
            alkalinity[i] = None
            emf0[i] = None
    dataset["alkalinity"] = pd.Series(alkalinity) * 1e6
    dataset["emf0"] = pd.Series(emf0)
    if verbose:
        print("Calkulate: solving complete!")
    return dataset
