import pandas as pd, numpy as np, PyCO2SYS as pyco2
from .. import convert, density, options, titrations


def get_measurement_type(dataset):
    if "measurement_type" not in dataset:
        dataset["measurement_type"] = "emf"
    dataset["measurement_type"] = np.where(
        np.isin(dataset.measurement_type, ["pH", "ph", "PH"]), "pH", "emf",
    )
    return dataset


def get_titrations(dataset, read_dat_kwargs=None):
    """(Re-)import all .dat files."""
    if "file_good" not in dataset:
        dataset["file_good"] = True
    if read_dat_kwargs is None:
        read_dat_kwargs = {}
    dataset.get_measurement_type()
    dats = {}
    for i, row in dataset.iterrows():
        if row.file_good:
            if "file_path" in row:
                fname = row.file_path + row.file_name
            else:
                fname = row.file_name
            try:
                dats[i] = titrations.read_dat(fname, **read_dat_kwargs)
                dats[i][row.measurement_type] = dats[i].measurement
            except IOError:
                print("Can't find file: '{}'.".format(fname))
                dats[i] = None
            except:
                print("Error importing file: '{}'.".format(fname))
                dats[i] = None
        else:
            dats[i] = None
    dataset["titration"] = pd.Series(dats)
    return dataset


def _get_analyte_temperature(row):
    if row.titration is not None:
        # Use the override temperature, if there is one
        if "temperature_override" in row:
            if ~np.isnan(row.temperature_override):
                analyte_temperature = row.temperature_override
            else:
                analyte_temperature = None
        else:
            analyte_temperature = None
        # If there isn't an override temperature, get it from the titration itself
        if analyte_temperature is None:
            try:
                analyte_temperature = float(row.titration.iloc[0].temperature)
            except ValueError:
                analyte_temperature = None
                print("Invalid temperature in file: '{}'".format(row.file_name))
    else:
        analyte_temperature = None
    return pd.Series({"analyte_temperature": analyte_temperature})


def get_analyte_temperature(dataset):
    dataset["analyte_temperature"] = dataset.apply(_get_analyte_temperature, axis=1)
    return dataset


def get_analyte_mass(dataset):
    """Calculate analyte mass in kg."""
    if "analyte_temperature" not in dataset:
        dataset.get_analyte_temperature()
    if "file_good" not in dataset:
        dataset["file_good"] = True
    dataset["file_good"] &= ~np.isnan(dataset.analyte_temperature)
    if "analyte_mass" not in dataset:
        dataset["analyte_mass"] = np.nan
        assert (
            "analyte_volume" in dataset
        ), "Dataset must contain either analyte_mass or analyte_volume."
    if "analyte_volume" in dataset:
        dataset["analyte_mass"] = np.where(
            np.isnan(dataset.analyte_mass),
            1e-3
            * dataset.analyte_volume
            * density.seawater_1atm_MP81(
                temperature=dataset.analyte_temperature, salinity=dataset.salinity
            ),
            dataset.analyte_mass,
        )
    return dataset


def get_titrant_density(dataset):
    if "molinity_HCl" not in dataset:
        dataset["molinity_HCl"] = options.molinity_HCl
    else:
        dataset["molinity_HCl"] = np.where(
            np.isnan(dataset.molinity_HCl), options.molinity_HCl, dataset.molinity_HCl
        )
    if "molinity_NaCl" not in dataset:
        dataset["molinity_NaCl"] = options.molinity_NaCl
    else:
        dataset["molinity_NaCl"] = np.where(
            np.isnan(dataset.molinity_NaCl),
            options.molinity_NaCl,
            dataset.molinity_NaCl,
        )
    if "titrant_density" not in dataset:
        dataset["titrant_density"] = np.nan
    dataset["titrant_density"] = np.where(
        np.isnan(dataset.titrant_density),
        density.HCl_NaCl_25C_DSC07(
            molinity_HCl=dataset.molinity_HCl, molinity_NaCl=dataset.molinity_NaCl
        ),
        dataset.titrant_density,
    )
    return dataset


def _get_titrant_mass(row):
    if row.titration is not None:
        try:
            if str(row.titrant_amount_unit).lower() == "g":
                row.titration["titrant_mass"] = row.titration.titrant_amount * 1e-3
            elif str(row.titrant_amount_unit).lower() == "kg":
                row.titration["titrant_mass"] = row.titration.titrant_amount
            elif str(row.titrant_amount_unit).lower() == "ml":
                row.titration["titrant_mass"] = (
                    row.titration.titrant_amount * row.titrant_density * 1e-3
                )
            else:
                row.titration["titrant_mass"] = (
                    row.titration.titrant_amount * row.titrant_density * 1e-3
                )
                print("Error: titrant_amount_unit must be either 'g' or 'ml'.")
        except TypeError:
            print("Invalid titration data in file: '{}'.".format(row.file_name))


def get_titrant_mass(dataset):
    """Calculate titrant mass in kg."""
    if "titrant_density" not in dataset:
        dataset.get_titrant_density()
    if "titrant_amount_unit" not in dataset:
        dataset["titrant_amount_unit"] = "ml"
    dataset.apply(_get_titrant_mass, axis=1)
    return dataset


def get_analyte_totals(dataset):
    """Make sure there is a molinity value for the undiluted analyte for every salt."""
    for salt in [
        "dic",
        "total_silicate",
        "total_phosphate",
        "total_ammonia",
        "total_sulfide",
        "total_alpha",
        "total_beta",
    ]:
        if salt not in dataset:
            dataset[salt] = 0
        dataset[salt] = np.where(np.isnan(dataset[salt]), 0, dataset[salt])
    if "opt_k_carbonic" not in dataset:
        dataset["opt_k_carbonic"] = options.pyco2_opt_k_carbonic
    else:
        dataset["opt_k_carbonic"] = np.where(
            np.isnan(dataset.opt_k_carbonic),
            options.pyco2_opt_k_carbonic,
            dataset.opt_k_carbonic,
        )
    if "opt_total_borate" not in dataset:
        dataset["opt_total_borate"] = options.pyco2_opt_total_borate
    else:
        dataset["opt_total_borate"] = np.where(
            np.isnan(dataset.opt_total_borate),
            options.pyco2_opt_total_borate,
            dataset.opt_total_borate,
        )
    totals = pyco2.salts.assemble(
        dataset.salinity,
        dataset.total_silicate,
        dataset.total_phosphate,
        dataset.total_ammonia,
        dataset.total_sulfide,
        dataset.opt_k_carbonic,
        dataset.opt_total_borate,
    )
    for salt in ["total_sulfate", "total_borate", "total_fluoride"]:
        if salt in dataset:
            dataset[salt] = np.where(
                np.isnan(dataset[salt]),
                totals[convert.calk_to_pyco2[salt]] * 1e6,
                dataset[salt],
            )
        else:
            dataset[salt] = totals[convert.calk_to_pyco2[salt]] * 1e6
    dataset["dic"] = np.where(np.isnan(dataset.dic), 0, dataset.dic)
    return dataset


def _get_titration_totals(row):
    if row.titration is not None:
        if "titrant_mass" in row.titration:
            row.titration["dilution_factor"] = convert.dilution_factor(
                row.analyte_mass, row.analyte_mass + row.titration.titrant_mass
            )
            for salt in [
                "dic",
                "total_silicate",
                "total_phosphate",
                "total_ammonia",
                "total_sulfide",
                "total_sulfate",
                "total_borate",
                "total_fluoride",
                "total_alpha",
                "total_beta",
            ]:
                row.titration[salt] = row[salt] * row.titration.dilution_factor


def get_titration_totals(dataset):
    """Get molinities of all salts throughout each titration."""
    dataset.apply(_get_titration_totals, axis=1)
    return dataset


def get_totals(dataset):
    """Get molinities of all salts in the analyte and throughout each titration."""
    dataset.get_analyte_totals()
    dataset.get_titration_totals()
    return dataset


def _get_k_constants(row):
    if row.titration is not None:
        if "titrant_mass" in row.titration:
            totals = {
                "Sal": row.salinity,
                "TNH3": row.titration.total_ammonia.values * 1e-6,
                "TPO4": row.titration.total_phosphate.values * 1e-6,
                "TSi": row.titration.total_silicate.values * 1e-6,
                "TH2S": row.titration.total_sulfide.values * 1e-6,
                "TB": row.titration.total_borate.values * 1e-6,
                "TF": row.titration.total_fluoride.values * 1e-6,
                "TSO4": row.titration.total_sulfate.values * 1e-6,
                "alpha": row.titration.total_alpha.values * 1e-6,
                "beta": row.titration.total_beta.values * 1e-6,
            }
            k_constants = {}
            for k in [
                "k_ammonia",
                "k_borate",
                "k_bisulfate",
                "k_carbonic_1",
                "k_carbonic_2",
                "k_fluoride",
                "k_phosphoric_1",
                "k_phosphoric_2",
                "k_phosphoric_3",
                "k_silicate",
                "k_sulfide",
                "k_water",
                "k_alpha",
                "k_beta",
            ]:
                if k in row:
                    if ~np.isnan(row[k]):
                        k_constants[convert.calk_to_pyco2[k]] = row[k]
            k_constants = pyco2.equilibria.assemble(
                row.titration.temperature,
                0,  # pressure
                totals,
                3,  # opt_pH_scale (i.e. Free)
                row.opt_k_carbonic,
                row.opt_k_bisulfate,
                row.opt_k_fluoride,
                row.opt_gas_constant,
                Ks=k_constants,
            )
            for k, v in k_constants.items():
                if k in convert.pyco2_to_calk:
                    row.titration[convert.pyco2_to_calk[k]] = v


def get_k_constants(dataset):
    if "analyte_temperature" not in dataset:
        dataset.get_analyte_temperature()
    if "opt_k_carbonic" not in dataset:
        dataset["opt_k_carbonic"] = options.pyco2_opt_k_carbonic
    else:
        dataset["opt_k_carbonic"] = np.where(
            np.isnan(dataset.opt_k_carbonic),
            options.pyco2_opt_k_carbonic,
            dataset.opt_k_carbonic,
        )
    if "opt_k_bisulfate" not in dataset:
        dataset["opt_k_bisulfate"] = options.pyco2_opt_k_bisulfate
    else:
        dataset["opt_k_bisulfate"] = np.where(
            np.isnan(dataset.opt_k_bisulfate),
            options.pyco2_opt_k_bisulfate,
            dataset.opt_k_bisulfate,
        )
    if "opt_k_fluoride" not in dataset:
        dataset["opt_k_fluoride"] = options.pyco2_opt_k_fluoride
    else:
        dataset["opt_k_fluoride"] = np.where(
            np.isnan(dataset.opt_k_fluoride),
            options.pyco2_opt_k_fluoride,
            dataset.opt_k_fluoride,
        )
    if "opt_gas_constant" not in dataset:
        dataset["opt_gas_constant"] = options.pyco2_opt_gas_constant
    else:
        dataset["opt_gas_constant"] = np.where(
            np.isnan(dataset.opt_gas_constant),
            options.pyco2_opt_gas_constant,
            dataset.opt_gas_constant,
        )
    for extra in ["alpha", "beta"]:
        k_extra = "k_" + extra
        if k_extra not in dataset:
            dataset[k_extra] = 1e-7
        else:
            dataset[k_extra] = np.where(
                np.isnan(dataset[k_extra]), 1e-7, dataset[k_extra],
            )
    dataset.apply(_get_k_constants, axis=1)
    return dataset
