# 2.1

Calkulate v2.1 is the first stable version in Python! New features will continue to be added, and bugs fixed, but the API of existing functions will no longer be changed.

## 2.1.0

**Release date:** forthcoming

  * Renamed many functions and variables, as now described in the [conventions](../conventions) documentation;
  * Added **plot** module to visualise titration results;
  * Updated `halfGran` solver to return values in dict field `'x'` for consistency;
  * Fixed fatal bug in `Dickson1981` solver;
  * Fixed incorrect `eqConstants['P1']` value in D81 dataset loading function;
  * Combined all calibration and VINDTA-style functions into generic functions that can implement any solving method;
  * Major documentation overhaul.

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
