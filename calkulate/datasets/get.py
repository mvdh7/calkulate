import pandas as pd
from .. import titrations


def get_titrations(dataset, **read_dat_kwargs):
    """(Re-)import all .dat files."""
    if "file_good" not in dataset:
        dataset["file_good"] = True
    dats = {}
    for i, row in dataset.iterrows():
        if row.file_good:
            if "file_path" in row:
                fname = row.file_path + row.file_name
            else:
                fname = row.file_name
            try:
                dats[i] = titrations.read_dat(fname, **read_dat_kwargs)
            except IOError:
                print("Can't find file: '{}'.".format(fname))
                dats[i] = None
            except:
                print("Error importing file: '{}'.".format(fname))
                dats[i] = None
    dataset["titration"] = pd.Series(dats)
    return dataset
