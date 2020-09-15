# Version history

## Version 3: Pythonic

Calkulate v3 switches to an object-oriented approach.  This makes working with individual titrations and handling large collections of them much less cumbersome.

[PyCO2SYS](https://PyCO2SYS.rtfd.io) is now used to determine equilibrium constants and estimate total salt concentrations from salinity, instead of replicating those functions in Calkulate.

### v3.0.0 (15 Sep 2020)

!!! info "Release notes for v3.0.0"

    * Introduces object-oriented approach in which most Calkulate functions are extensions to pandas DataFrames.
    * Uses PyCO2SYS to evaluate all equilibrium constants and salt concentrations from salinity and temperature.

## Version 2: Python

Calkulate v2 transitioned from MATLAB to Python and added several additional alkalinity solvers.  Despite being implemented in Python, the style was still rather MATLAB-esque.

## Version 1: MATLAB

The original Calkulate v1 was a MATLAB implementation of the "half-Gran" alkalinity solver described by [H15](../references/#h). The final version (v1.0.2) remains [freely available](https://github.com/mvdh7/calkulate/tree/1.0.2).
