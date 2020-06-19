# Input data formatting

To work with a dataset of multiple titrations you need to provide two things:

  1.  A [**titration table**](#the-titration-table): the metadata for each titration.
  2.  The [**titration files**](#titration-files): a text file for each titration containing the measurements. 

## The titration table

The titration table is a [pandas DataFrame](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html) or anything that can be converted into one, such as a dict, CSV file, or Excel spreadsheet.

Each row contains the metadata for a single titration.  A very simple CSV file for a single titration could look like this:

    file_name          , salinity , analyte_volume
    titration_file.dat , 35       , 100

### Titration table columns

Click on each banner to expand and see more detail.

Extra columns will be ignored by Calkulate and should not cause any problems.

#### Required columns

??? success "`file_name`: *name of the titration file*"
    This should include the file extension (e.g. `.dat`).  The file path can either be included in the file name or be specified separately with `file_path`.

??? success "`salinity`: *practical salinity of the analyte*"
    This is dimensionless.

??? question "`analyte_mass` or `analyte_volume`: *mass or volume of the analyte*"
    Mass must be in g and volume in ml.
    
    If volume is provided, it is converted into mass following [MP81](../references/#m), assuming that the acid contains 0.6 M NaCl and 0.1 M HCl.
    
    If both are provided for a given titration, mass is used and volume is ignored.

#### Optional columns

If optional column values are only needed for some of your titrations, just use `np.nan` in rows where they are not required.

??? tip "`alkalinity_certified`: *known alkalinity value (e.g. for a reference material)*"
    Must be in μmol/kg-solution.

??? tip "`analysis_batch`: *identifies subsets that can be calibrated together*"
    A list of analysis batch names or numbers.  It is assumed that all titrations with the same `analysis_batch` value have a common titrant.  If nothing is provided, all titrations are assumed to have a common titrant.

??? tip "`bisulfate_constant`: *which bisulfate dissociation constant to use*"
    Sets whether the bisulfate dissociation constant is estimated from salinity and temperature by [PyCO2SYS](https://PyCO2SYS.rtfd.io) following:
    
    * `1`: [D90a](https://PyCO2SYS.readthedocs.io/en/latest/refs/#d) (default), or
    * `2`: [KRCB77](https://PyCO2SYS.readthedocs.io/en/latest/refs/#k).

??? tip "`borate_ratio`: *which total borate:salinity relationship to use*"
    Sets whether total borate is estimated from salinity by [PyCO2SYS](https://PyCO2SYS.rtfd.io) following:
    
    * `1`: [U74](https://PyCO2SYS.readthedocs.io/en/latest/refs/#u), or
    * `2`: [LKB10](https://PyCO2SYS.readthedocs.io/en/latest/refs/#l) (default).

??? tip "`carbonic_constants`: *which carbonic acid dissociation constants to use*"
    Sets which carbonic acid constants to use from PyCO2SYS, can be any integer from `1` to `15` inclusive.  Default is `10` for [LDK00](https://pyco2sys.readthedocs.io/en/latest/refs/#l).  See the [PyCO2SYS docs on `K1K2CONSTANTS`](https://pyco2sys.readthedocs.io/en/latest/co2sys/#settings) for details.

??? tip "`file_good`: *is the titration file valid?*"
    Where set to `False`, Calkulate does not attempt to import the corresponding titration file.

??? tip "`file_path`: *path to the titration file*"
    This is simply prefixed to `file_name`.

??? tip "`fluoride_constant`: *which HF dissociation constant to use*"
    Sets whether the HF dissociation constant is estimated from salinity and temperature by [PyCO2SYS](https://PyCO2SYS.rtfd.io) following:
    
    * `1`: [DR79](https://PyCO2SYS.readthedocs.io/en/latest/refs/#d) (default), or
    * `2`: [PF87](https://PyCO2SYS.readthedocs.io/en/latest/refs/#p).

??? tip "`k_*`: *stoichiometric equilibrium constants*"
    If not provided, then these are calculated from temperature and salinity by [PyCO2SYS](https://PyCO2SYS.rtfd.io).

    * `k_boric`: boric acid equilibrium constant.
    * `k_carbonic_1` and `k_carbonic_2`: carbonic acid dissociation constants.
    * `k_hydrofluoric`: HF equilibrium constant.
    * `k_phosphoric_1`, `k_phosphoric_2` and `k_phosphoric_3`: phosphoric acic dissociation constants.
    * `k_orthosilicic`: silicate equilibrium constant.
    * `k_hydrosulfuric_1`: hydrogen sulfide equilibrium constant.
    * `k_sulfuric_2`: bisulfate dissociation constant.
    * `k_water`: water equilibrium constant.

??? tip "`measurement_type`: *type of measurement in the titration file*"
    Use `"EMF"` (default) for potentiometric measurements, or `"pH"` for direct pH measurements.

??? tip "`reference_good`": *use this reference material to calibrate?*
    Set to `False` for reference materials that you do not wish to include in determining batch-averaged titrant concentrations.

??? tip "`temperature_override`: *titration temperature*"
    Must be in °C.  If not supplied, the values in the titration file are used.  Otherwise, this value overrides them.

??? tip "`titrant_amount_unit`: *unit for the amount of titrant in the titration file*"
    Either `"ml"` (default) for volume in ml or `"g"` for mass in g.  If volume is provided, it is converted into mass following [DSC07](../references/#d), assuming the acid contains 0.6 M NaCl and 0.1 M HCl.

??? tip "`titrant_concentration`: *concentration of the titrant*"
    Must be in mol/l.  Is converted to mol/kg-solution following [DSC07](../references/#d).

??? tip "`titrant_molinity`: *molinity of the titrant*"
    Must be in mol/kg-solution.

??? tip "`total_*`: *total salt concentrations*"
    These must all be in μmol/kg-solution.

    Some default to zero if not provided:

    * `total_ammonia`: total ammonia + ammonium
    * `total_carbonate`: total dissolved inorganic carbon
    * `total_phosphate`: total phosphate
    * `total_silicate`: total silicate
    * `total_sulfide`: total hydrogen sulfide

    Others are estimated from salinity by [PyCO2SYS](https://PyCO2SYS.rtfd.io) if not provided:

    * `total_borate`: total borate
    * `total_fluoride`: total fluoride
    * `total_sulfate`: total sulfate

## Titration files

Each titration file is a text file containing the measurements of the solution carried out during a titration.  The file should contain data in columns, where each row represents a measurement after a separate titrant addition.  There must be at least three columns, containing:

  1. The amount of titrant added to the analyte in ml or g.
  2. The EMF measured across the titrant-analyte mixture in mV, or its pH.
  3. The temperature of the titrant-analyte mixture in °C.

By default, Calkulate expects that:

  * There are two lines at the start of the file to be ignored before the columns begin.
  * The above are the first three columns to appear in the text file, and in the order given.
  * The file is tab-delimited.
  * Nothing comes after the columns of titration data in the file.

For example, a file could contain the following:

    This first line is ignored.
    This second line is also ignored.
    0.00    183.1   25.2
    0.50    225.4   25.1
    1.00    290.3   25.0
    1.50    343.4   25.1

However, all these assumptions can be adjusted, as described in [Individual titrations](../titrations).