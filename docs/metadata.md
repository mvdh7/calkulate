# The Dataset of titration metadata

A Calkulate Dataset is just a [pandas DataFrame](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html) with [some extra methods added](methods).

You can make a dataset [from an existing DataFrame](../io/#convert-from-a-dataframe) or by [importing a spreadsheet from various formats](../io/#import-from-excel-csv-or-dbs).

Whichever way you do it, each row of the Dataset corresponds to a separate titration, and each column contains a different piece of metadata.  Some columns are required, others are optional.  The columns must be named as follows.

## Dataset column names

Click on each banner to expand and see more detail.

Extra columns will be ignored by Calkulate and should not cause any problems.

### Required columns

??? success "`file_name` : *name of the titration file*"
    This should include the file extension (e.g. `.dat`).  The file path can either be included in the file name or be specified separately in a `file_path` column.

??? success "`salinity` : *practical salinity of the analyte*"
    This is dimensionless.

??? question "`analyte_mass` or `analyte_volume` : *mass or volume of the analyte*"
    **Mass must be in kg and volume in ml.**
    
    If volume is provided, it is converted into mass following [MP81](../references/#m).
    
    If both are provided for a given titration, mass is used and volume is ignored.

### Recommended columns

Not strictly essential, but you'll need these most of the time.  If values are only needed for some of your titrations, just use `np.nan` in rows where they are not required.

??? tip "`alkalinity_certified` : *known alkalinity value (e.g. for a reference material)*"
    Must be in μmol/kg-solution.

??? tip "`analysis_batch` : *identifies subsets that can be calibrated together*"
    A list of analysis batch names or numbers.  It is assumed that all titrations with the same `analysis_batch` value have a common titrant.  If nothing is provided, all titrations are assumed to have a common titrant.

??? tip "`dic` : *dissolved inorganic carbon*"
    In μmol/kg-solution.

    Defaults to zero if not provided.

??? tip "`file_path` : *path to the titration file*"
    This is prefixed to `file_name`.

??? tip "`total_*` : *total salt concentrations*"
    These must all be in μmol/kg-solution.

    Some default to zero if not provided:

    * `total_ammonia` : total ammonia + ammonium
    * `total_phosphate` : total phosphate
    * `total_silicate` : total silicate
    * `total_sulfide` : total hydrogen sulfide

    Others are estimated from salinity by [PyCO2SYS](https://PyCO2SYS.rtfd.io) if not provided:

    * `total_borate` : total borate
    * `total_fluoride` : total fluoride
    * `total_sulfate` : total sulfate

### Optional columns

Like for the recommended columns, if optional column values are only needed for some of your titrations, just use `np.nan` in rows where they are not required.

??? info "`bisulfate_constant` : *which bisulfate dissociation constant to use*"
    Sets whether the bisulfate dissociation constant is estimated from salinity and temperature by [PyCO2SYS](https://PyCO2SYS.rtfd.io) following:
    
    * `1` : [D90a](https://PyCO2SYS.readthedocs.io/en/latest/refs/#d) (default), or
    * `2` : [KRCB77](https://PyCO2SYS.readthedocs.io/en/latest/refs/#k).

??? info "`borate_ratio` : *which total borate:salinity relationship to use*"
    Sets whether total borate is estimated from salinity by [PyCO2SYS](https://PyCO2SYS.rtfd.io) following:
    
    * `1` : [U74](https://PyCO2SYS.readthedocs.io/en/latest/refs/#u), or
    * `2` : [LKB10](https://PyCO2SYS.readthedocs.io/en/latest/refs/#l) (default).

??? info "`carbonic_constants` : *which carbonic acid dissociation constants to use*"
    Sets which carbonic acid constants to use from PyCO2SYS, can be any integer from `1` to `15` inclusive.  Default is `10` for [LDK00](https://pyco2sys.readthedocs.io/en/latest/refs/#l).  See the [PyCO2SYS docs on `K1K2CONSTANTS`](https://pyco2sys.readthedocs.io/en/latest/co2sys/#settings) for details.

??? info "`file_good` : *is the titration file valid?*"
    Where set to `False`, Calkulate does not attempt to import the corresponding titration file.

??? info "`fluoride_constant` : *which HF dissociation constant to use*"
    Sets whether the HF dissociation constant is estimated from salinity and temperature by [PyCO2SYS](https://PyCO2SYS.rtfd.io) following:
    
    * `1` : [DR79](https://PyCO2SYS.readthedocs.io/en/latest/refs/#d) (default), or
    * `2` : [PF87](https://PyCO2SYS.readthedocs.io/en/latest/refs/#p).

??? info "`k_*` : *stoichiometric equilibrium constants*"
    If not provided, then these are calculated from temperature and salinity by [PyCO2SYS](https://PyCO2SYS.rtfd.io).

    They must all be on the Free pH scale.

    * `k_borate` : boric acid equilibrium constant.
    * `k_carbonic_1` and `k_carbonic_2` : carbonic acid dissociation constants.
    * `k_fluoride` : HF equilibrium constant.
    * `k_phosphoric_1`, `k_phosphoric_2` and `k_phosphoric_3` : phosphoric acic dissociation constants.
    * `k_orthosilicic` : silicate equilibrium constant.
    * `k_sulfide` : hydrogen sulfide equilibrium constant.
    * `k_bisulfate` : bisulfate dissociation constant.
    * `k_water` : water equilibrium constant.

??? info "`measurement_type` : *type of measurement in the titration file*"
    Use `"emf"` (default) for potentiometric measurements, or `"pH"` for direct pH measurements.

??? info "`molinity_HCl` : *approximate HCl molinity in the titrant*"
    In mol/kg-solution.  Defaults to 0.1 mol/kg-solution if not provided.  Only used to estimate titrant density, not for calibration.

??? info "`molinity_NaCl` : *approximate NaCl molinity in the titrant*"
    In mol/kg-solution.  Defaults to 0.6 mol/kg-solution if not provided.  Used to estimate titrant density.

??? info "`reference_good` : *use this reference material to calibrate?*"
    Set to `False` for reference materials that you do not wish to include in determining batch-mean titrant concentrations.

??? info "`temperature_override` : *titration temperature*"
    Must be in °C.  If not supplied, the values in the titration file are used.  Otherwise, this value overrides them.

??? info "`titrant_amount_unit` : *unit for the amount of titrant in the titration file*"
    Either `"ml"` (default) for volume in ml or `"g"` / `"kg"` for mass in g / kg.
    
    If volume is provided, it is converted into mass following [DSC07](../references/#d), assuming the acid contains 0.6&nbsp;M NaCl and 0.1&nbsp;M HCl.  If you used a different ratio, use the columns `molinity_NaCl` and `molinity_HCl` to set this.

??? info "`titrant_concentration` : *concentration of the titrant*"
    Must be in mol/l.  Is converted to mol/kg-solution following [DSC07](../references/#d).

??? info "`titrant_molinity` : *molinity of the titrant*"
    Must be in mol/kg-solution.
