# Version history

## Version 23: best of v2 and v3

Calkulate v3 went too far overboard with the OO approach and ended up being very slow and overly complex behind the scenes as a result.  Calkulate v23 therefore mashes together the best bits of v2 and v3 for the ultimate alkalinity solving experience.

### 23.6 (19 February 2024)

!!! info "Changes in v23.6"

    * Added missing components to titration table in a `Titration`.

    ***v23.6.1 bug fixes (20 April 2024)***

    * Removed excessive `print` statements from debugging.


### 23.5 (4 July 2023)

!!! info "Changes in v23.5"

    * Added `to_pandas` method to `Dataset` for a more convenient downgrade to a pandas `DataFrame`.
    * Improved docstrings for several key functions and classes.
    * Fixed bugs in logicals in various `calk.dataset` functions which had meant that NaN values in the metadata file were not always being ignored.
    * Fixed bug in the `Dataset.solve` method such that the `analyte_mass` is now automatically calculated from `analyte_volume`.
    * Switched to using pyproject.toml instead of setup.py for building the package.
    * Renamed module `io` as `read`.

### 23.4 (13 June 2023)

!!! info "Changes in v23.4"

    * Swapped deprecated pandas `.iteritems` for `.items`.
    * Improvements to `Titration` class.

### 23.3 (22 June 2022)

!!! info "Changes in v23.3"

    * Added a test of the [DSC07](../references/#d) SOP 3b example dataset.
    * Updated for compatibility with PyCO2SYS v1.8.1 (but no longer with v1.8.0).

### 23.2

Adds `Titration` class for investigating single titrations.

#### 23.2.2 (4 February 2022)

!!! info "Changes in v23.2.2"

    * More properties calculated during titrations.

#### 23.2.1 (13 August 2021)

!!! info "Changes in v23.2.1"

    * Solver functions return additional diagnostic information about each titration.

#### 23.2.0 (6 July 2021)

!!! info "Changes in v23.2.0"

    * Adds (as yet undocumented) `Titration` class for investigating single titrations.
    * Adds (as yet undocumented) plotting functions for `Titration` objects.

### 23.1

Adds solvers for titrations with an H<sub>2</sub>SO<sub>4</sub> titrant.

#### 23.1.1 (7 June 2021)

!!! info "Changes in v23.1.1"

   * Now compatible with PyCO2SYS v1.7.0.

#### 23.1.0 (29 March 2021)

!!! info "Changes in v23.1.0"

    * Adds solvers for titrations with an H<sub>2</sub>SO<sub>4</sub> titrant (see optional `titrant` column in the [metadata contents](../metadata/#optional-columns)).
    * Adds optional `titrant_density` column to overwrite the internally calculated titrant density with a user-specified value.

### 23.0

#### 23.0.2 (15 March 2021)

!!! info "Changes in v23.0.2"

    * Better handling of missing dates in files imported with `read_dbs`.
    * Minor adjustments and bug fixes in a few internal functions.

#### 23.0.1 (25 February 2021)

!!! info "Changes in v23.0.1"

    * Print more informative error messages when titration data files cannot be found.
    * Use `read_dat_genfromtxt` by default if `read_dat_method` not recognised, instead of throwing an error.

#### 23.0.0 (22 February 2021)

!!! info "Introduction of v23.0.0"

    * Object-oriented syntax from v3 is available to quickly work with datasets of many titrations at once.
    * Underlying functions work much faster with raw NumPy arrays as in v2.
    * Some basic plotting functions added.

## Version 3: poorly Pythonic

Calkulate v3 switches to an object-oriented approach.  This makes working with individual titrations and handling large collections of them much less cumbersome.

[PyCO2SYS](https://PyCO2SYS.rtfd.io) is now used to determine equilibrium constants and estimate total salt concentrations from salinity, instead of replicating those functions in Calkulate.

### 3.1

#### 3.1.0 (27 October 2020)

!!! info "Changes in v3.1.0"

    ***Major bug fix***

    * Fixed unit conversion bug when evaluating equilibrium constants in PyCO2SYS.  **All results calculated using v3.0.X should be redetermined!**

    ***Better consistency with PyCO2SYS***

    * Updated for compatability with PyCO2SYS v1.6.0.
    * Added two optional extra alkalinity components.
    * Renamed various internal variables for better consistency.

### 3.0

#### 3.0.1 (23 September 2020)

!!! info "Changes in v3.0.1"

    ***Bug fixes***

    * Fix bug in handling `data["file_good"] = False` cases.
    * Set DIC to zero where its value is NaN.
    * Skip over errors in titration files with a warning rather than throw a breaking error.

#### 3.0.0 (15 September 2020)

!!! info "Release notes for v3.0.0"

    * Introduces object-oriented approach in which most Calkulate functions are extensions to pandas DataFrames.
    * Uses PyCO2SYS to evaluate all equilibrium constants and salt concentrations from salinity and temperature.

## Version 2: Python

Calkulate v2 transitioned from MATLAB to Python and added several additional alkalinity solvers.  Despite being implemented in Python, the style was still rather MATLAB-esque.

## Version 1: MATLAB

The original Calkulate v1 was a MATLAB implementation of the "half-Gran" alkalinity solver described by [H15](../references/#h). The final version (v1.0.2) remains [freely available](https://github.com/mvdh7/calkulate/tree/1.0.2).
