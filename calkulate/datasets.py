import pandas as pd, numpy as np
from . import density, io


def has_value(row, col):
    """Check if a row exists and if it is not NaN."""
    has_value = col in row
    if has_value:
        has_value &= ~np.isnan(row[col])
    return has_value


def get_dat_files(dataset, **read_dat_kwargs):
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
                dats[i] = io.read_dat(fname, **read_dat_kwargs)
            except IOError:
                print("Can't find file: '{}'.".format(fname))
                dats[i] = None
            except:
                print("Error importing file: '{}'.".format(fname))
                dats[i] = None
    dataset["dat_dict"] = pd.Series(dats)
    return dataset


def prepare_row(
    dataset,
    ix,
    titrant_density_func=density.HCl_NaCl_25C_DSC07,
    titrant_density_kwargs=None,
):
    """Prepare Calkulate inputs from a Dataset row."""
    row = dataset.loc[ix]
    titrant = {}
    analyte = {}
    temperature = row.dat_dict["mixture_temperature"]
    if has_value(row, "temperature_override"):
        temperature[:] = row.temperature_override
    titrant["density"] = titrant_density_func(**titrant_density_kwargs)
    if "titrant_amount_unit" in row:
        if row.titrant_amount_unit == "g":
            amount_factor = 1
        elif row.titrant_amount_unit == "ml":
            amount_factor = 1e-3 / titrant["density"]
    else:  # assume ml by default
        amount_factor = 1e-3 / titrant["density"]
    titrant["mass"] = row.dat_dict["titrant_amount"] * amount_factor
    if has_value(row, "titrant_concentration"):
        titrant["molinity"] = row.titrant_molinity * 1e-3 / titrant["density"]
    if has_value(row, "titrant_molinity"):
        titrant["molinity"] = row.titrant_molinity
    if has_value(row, "alkalinity_certified"):
        analyte["alkalinity_certified"] = row.alkalinity_certified
    if has_value(row, "dic"):
        analyte["dic"] = row.dic
    else:
        analyte["dic"] = 0
    if has_value(row, "analyte_mass"):
        analyte["mass"] = row.analyte_mass
    elif has_value(row, "analyte_volume"):
        analyte["mass"] = (
            row.analyte_volume
            * 1e-3
            / density.seawater_1atm_MP81(
                temperature=temperature[0], salinity=row.salinity
            )
        )
    else:
        print("Error: you need to provide either analyte_mass or analyte_volume!")
    return titrant, analyte, emf_or_pH, temperature, measurement_type


class Dataset(pd.DataFrame):
    get_dat_files = get_dat_files
    prepare_row = prepare_row
