# Input data formatting

To work with titration data in Calkulate you need to provide two things:

  1.  The [titration file(s)](#titration-files): a text file for each titration containing the measurements. 
  2.  A [titration table](#the-titration-table): the metadata for each titration.

Even if you are only using data from a single titration, you still need to make a single-row titration table to work with it in Calkulate.

## Titration files

Each titration file is a text file containing the measurements of the solution carried out during a titration.  The file should contain data in columns, where each row represents a measurement after a separate titrant addition.  There must be at least three columns:

  1. The amount of titrant added to the analyte in ml or g.
  2. The EMF measured across the titrant-analyte mixture in mV, or its pH.
  3. The temperature of the titrant-analyte mixture in °C.

By default, Calkulate expects that:

  * There are two lines at the start of the file to be ignored before the data columns described above begin.
  * The columns appear next, in the order given (from left to right).
  * The columns are tab-delimited.
  * Nothing comes after the columns of titration data in the file.

For example, a file could contain the following:

    This first line is ignored.
    This second line is also ignored.
    0.00    183.1   25.2
    0.50    225.4   25.1
    1.00    290.3   25.0
    1.50    343.4   25.1

!!! tip "Titration files in different formats"

    #### Import settings

    All assumptions above about what your titration files look like can be adjusted.  The first thing you should do is work out what adjustments you need to make, if any.  You should repeat this each time you have a titration file in a new format.

    If your titration data files come from a VINDTA, you should be able to skip this step.

    Internally, Calkulate imports titration files using:

        :::python
        import calkulate as calk
        tt = calk.io.read_dat(
            fname,
            titrant_amount_col=0,
            measurement_col=1,
            temperature_col=2,
            delimiter="\t",
            skip_header=2,
            **kwargs
        )

    The only required input, `fname`, is the titration file name.  Calkulate imports the data from this file using [`numpy.genfromtxt`](https://numpy.org/doc/stable/reference/generated/numpy.genfromtxt.html), which is where inputs `delimiter`, `skip_header`, and any other `kwargs` that you may need to use are passed on to.

    `np.genfromtxt` returns all the data in the titration file as a 2-dimensional NumPy array.  `calkulate.io.read_dat` then extracts the columns identified by `titrant_amount_col`, `measurement_col` and `temperature_col` into fields `"titrant_amount"`, `"mixture_measurement"` and `"mixture_temperature"` of a dict `tt`.

    Before going further with Calkulate, you should make sure that you can import one of your titration files successfully using `calk.io.read_dat`.  You won't actually call this function directly, but you need to know any optional settings it requires for your titration files.  Assemble these optional settings in a dict and check they work.  For example, if you needed to skip 5 header lines, not 2, your file was comma-delimited, not tab-, and there was an extra column between the measurement data and the temperature:

    ```python
    import calkulate as calk
    
    # Assemble input arguments
    fname = "my_dat_file.dat"
    read_dat_kwargs = {
        "skip_header": 5,
        "temperature_col": 3,
        "delimiter": ",",
    }

    # Run the import function
    tt = calk.io.read_dat(fname, **read_dat_kwargs)
    
    # Check the output looks good
    print(tt["titrant_amount")  # prints out the amount of titrant column
    print(tt["mixture_measurement"])  # prints out the measurement column (e.g. EMF)
    print(tt["mixture_temperature"])  # prints out the temperature
    ```

    Keep a note of `read_dat_kwargs` - you'll need it later!

## The titration table

The titration table is a [pandas DataFrame](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html) or anything that can be converted into one, such as a dict, CSV file, or Excel spreadsheet.

Each row contains the metadata for a single titration.  A very simple CSV file for a single titration could look like this:

    file_name          , salinity , analyte_volume
    titration_file.dat , 35       , 100

Alternatively, you could define the same information in a dict:

    :::python
    titration_table = {
        "file_name": "titration_file.dat",
        "salinity": 35,
        "analyte_volume": 100,
    }

### Titration table columns

Click on each banner to expand and see more detail.

Extra columns will be ignored by Calkulate and should not cause any problems.

#### Required columns

??? success "`file_name`: *name of the titration file*"
    This should include the file extension (e.g. `.dat`).  The file path can either be included in the file name or be specified separately in a `file_path` column.

??? success "`salinity`: *practical salinity of the analyte*"
    This is dimensionless.

??? question "`analyte_mass` or `analyte_volume`: *mass or volume of the analyte*"
    Mass must be in g and volume in ml.
    
    If volume is provided, it is converted into mass following [MP81](../references/#m).
    
    If both are provided for a given titration, mass is used and volume is ignored.

#### Recommended columns

Not strictly essential, but you'll need these most of the time.  If values are only needed for some of your titrations, just use `np.nan` in rows where they are not required.

??? tip "`alkalinity_certified`: *known alkalinity value (e.g. for a reference material)*"
    Must be in μmol/kg-solution.

??? tip "`analysis_batch`: *identifies subsets that can be calibrated together*"
    A list of analysis batch names or numbers.  It is assumed that all titrations with the same `analysis_batch` value have a common titrant.  If nothing is provided, all titrations are assumed to have a common titrant.

??? tip "`file_path`: *path to the titration file*"
    This is simply prefixed to `file_name`.

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

#### Optional columns

Like for the recommended columns, if optional column values are only needed for some of your titrations, just use `np.nan` in rows where they are not required.

??? info "`bisulfate_constant`: *which bisulfate dissociation constant to use*"
    Sets whether the bisulfate dissociation constant is estimated from salinity and temperature by [PyCO2SYS](https://PyCO2SYS.rtfd.io) following:
    
    * `1`: [D90a](https://PyCO2SYS.readthedocs.io/en/latest/refs/#d) (default), or
    * `2`: [KRCB77](https://PyCO2SYS.readthedocs.io/en/latest/refs/#k).

??? info "`borate_ratio`: *which total borate:salinity relationship to use*"
    Sets whether total borate is estimated from salinity by [PyCO2SYS](https://PyCO2SYS.rtfd.io) following:
    
    * `1`: [U74](https://PyCO2SYS.readthedocs.io/en/latest/refs/#u), or
    * `2`: [LKB10](https://PyCO2SYS.readthedocs.io/en/latest/refs/#l) (default).

??? info "`carbonic_constants`: *which carbonic acid dissociation constants to use*"
    Sets which carbonic acid constants to use from PyCO2SYS, can be any integer from `1` to `15` inclusive.  Default is `10` for [LDK00](https://pyco2sys.readthedocs.io/en/latest/refs/#l).  See the [PyCO2SYS docs on `K1K2CONSTANTS`](https://pyco2sys.readthedocs.io/en/latest/co2sys/#settings) for details.

??? info "`file_good`: *is the titration file valid?*"
    Where set to `False`, Calkulate does not attempt to import the corresponding titration file.

??? info "`fluoride_constant`: *which HF dissociation constant to use*"
    Sets whether the HF dissociation constant is estimated from salinity and temperature by [PyCO2SYS](https://PyCO2SYS.rtfd.io) following:
    
    * `1`: [DR79](https://PyCO2SYS.readthedocs.io/en/latest/refs/#d) (default), or
    * `2`: [PF87](https://PyCO2SYS.readthedocs.io/en/latest/refs/#p).

??? info "`k_*`: *stoichiometric equilibrium constants*"
    If not provided, then these are calculated from temperature and salinity by [PyCO2SYS](https://PyCO2SYS.rtfd.io).

    * `k_boric`: boric acid equilibrium constant.
    * `k_carbonic_1` and `k_carbonic_2`: carbonic acid dissociation constants.
    * `k_hydrofluoric`: HF equilibrium constant.
    * `k_phosphoric_1`, `k_phosphoric_2` and `k_phosphoric_3`: phosphoric acic dissociation constants.
    * `k_orthosilicic`: silicate equilibrium constant.
    * `k_hydrosulfuric_1`: hydrogen sulfide equilibrium constant.
    * `k_sulfuric_2`: bisulfate dissociation constant.
    * `k_water`: water equilibrium constant.

??? info "`measurement_type`: *type of measurement in the titration file*"
    Use `"EMF"` (default) for potentiometric measurements, or `"pH"` for direct pH measurements.

??? info "`reference_good`: *use this reference material to calibrate?*"
    Set to `False` for reference materials that you do not wish to include in determining batch-averaged titrant concentrations.

??? info "`temperature_override`: *titration temperature*"
    Must be in °C.  If not supplied, the values in the titration file are used.  Otherwise, this value overrides them.

??? info "`titrant_amount_unit`: *unit for the amount of titrant in the titration file*"
    Either `"ml"` (default) for volume in ml or `"g"` for mass in g.  If volume is provided, it is converted into mass following [DSC07](../references/#d), assuming the acid contains 0.6 M NaCl and 0.1 M HCl.

??? info "`titrant_concentration`: *concentration of the titrant*"
    Must be in mol/l.  Is converted to mol/kg-solution following [DSC07](../references/#d).

??? info "`titrant_molinity`: *molinity of the titrant*"
    Must be in mol/kg-solution.
