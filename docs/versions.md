# 2.3

Calkulate v2.3 introduces a new `Potentiometric` class that titration datasets can be imported into, enabling easier high-level analysis using the existing array-based functions. This is not yet properly documented, although there is a working example script in the Github repo (*examples/stepwise-alkalinity.py*).

## 2.3.0

**Release date:** 2020-04-02

  * Added `convert` module with convenience functions to convert acid and sample volumes into masses.
  * Relocated EMF to [H<sup>+</sup>] conversion functions from `solve` to `convert`, with aliases to maintain backwards compatibility.
  * Added `titration` module with new `Potentiometric` class for higher-level manipulation of titration datasets.
  * Added plotting functions for `Potentiometric` titration objects into the `plot` module.

---

# 2.2

The main advance in Calkulate v2.2 is switching to using [PyCO2SYS](https://github.com/mvdh7/PyCO2SYS) to evaluate all equilibrium constants and solute concentrations from temperature and salinity, instead of having an independent set of equivalent functions written into Calkulate itself.

**Release date:** 2020-03-20

## 2.2.0

  * Updated `concentrations.concTotals` and `dissociation.eqConstants` to use [PyCO2SYS](https://github.com/mvdh7/PyCO2SYS) functions to calculate things instead of functions built in to Calkulate.
  * Removed all equilibrium constant and concentration functions that were previously internally in Calkulate.
  * Added optional total ammonia and hydrogen sulfide concentration inputs. If values are provided, these equilibria will now be taken into account by the `complete` solver.
  * Fixed `io.writeDat` function to correctly delimit .dat files with tabs (not spaces).
  * Added `simulate.titration` function to directly simulate a titration dataset.
  * Updated solver functions to optionally accept `concTotals['C']` as an array.
  * Renamed module `vindta` to `datfile`, but added alias to avoid breaking existing code.

---

# 2.1

Calkulate v2.1 is the first stable version in Python! New features will continue to be added, documentation developed, and bugs fixed, but the API of existing functions will no longer be changed.

## 2.1.0

**Release date:** 2019-08-06

  * Renamed many functions and variables, as now described in the [conventions](../conventions) documentation;
  * Added **plot** module to visualise titration results;
  * Updated `halfGran` solver to return values in dict field `'x'` for consistency;
  * Fixed fatal bug in `Dickson1981` solver;
  * Fixed incorrect `eqConstants['P1']` value in D81 dataset loading function;
  * Combined all calibration and VINDTA-style functions into generic functions that can implement any solving method;
  * Major documentation overhaul, and added workflow examples.

---

# 2.0 (Python, beta)

Despite the version number, Calkulate v2.0 is still a beta, as the original MATLAB program is extended and converted to Python.

## 2.0.22

**Release date:** 2019-04-09

  * Added full outline documentation and function docstrings;
  * Added **meta** module, with single-source-of-truth version number;
  * Converted lists of solute concentrations (`XT`) and dissociation constants (`KX`) to dicts;
  * Renamed *MPH* functions to *complete*;
  * Updated various aspects of nomenclature throughout;
  * Switched pH simulator to solve for pH, not [H<sup>+</sup>];
  * Active efforts towards MATLAB integration paused.

## 2.0.21

**Release date:** 2019-01-15

  * All functionality from original MATLAB version converted to Python, with equivalent results;
  * Tested all solvers successfully against simulated titration data from Dickson (1981);
  * Final release before starting to track updates here in the version history.

---

# 1.0.2 (MATLAB)

**Release date:** 2015

The original Calkulate v1.0.2 is a MATLAB-only implementation of the half-Gran alkalinity solver, as described by Humphreys (2015). The final version will remain freely available: see [the appropriate project branch on GitHub](https://github.com/mvdh7/calkulate/tree/1.0.2).
