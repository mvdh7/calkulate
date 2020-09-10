import pandas as pd
from . import io


def get_dat_files(dbs, **read_dat_kwargs):
    """(Re-)import all .dat files."""
    if "file_good" not in dbs:
        dbs["file_good"] = True
    dats = {}
    for i, row in dbs.iterrows():
        if row.file_good:
            if "file_path" in row:
                fname = row.file_path + row.file_name
            else:
                fname = row.file_name
            try:
                dats[i] = io.read_dat(fname, **read_dat_kwargs)
            except IOError:
                print("Can't find file: '{}'.".format(fname))
                dats[i] = None
            except:
                print("Error importing file: '{}'.".format(fname))
                dats[i] = None
    dbs["dat_dict"] = pd.Series(dats)
    return dbs


class Dataset(pd.DataFrame):
    get_dat_files = get_dat_files
