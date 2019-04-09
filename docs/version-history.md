# 2.0 (Python, beta)

Despite the version number, Calkulate v2.0 is in beta, as the original MATLAB program is being extended and converted to Python.


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

<hr />

# 1.0.2 (MATLAB)

**Release date:** 2015

The original Calkulate v1.0.2 is a MATLAB-only implementation of the half-Gran alkalinity solver, as described by Humphreys (2015). The final version will remain freely available: see [the appropriate branch on GitHub](https://github.com/mvdh7/calkulate/tree/1.0.2).
