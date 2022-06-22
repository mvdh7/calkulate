# Validating Calkulate

## Automatic tests

Several of the [automatic tests](https://github.com/mvdh7/calkulate/tree/main/tests) include checks on whether Calkulate does self-consistent calibrations, for example solving a sample where the titrant was calibrated with itself, which should return exactly the correct alkalinity.

So if the Tests badge below says "passing":

![Tests](https://github.com/mvdh7/calkulate/workflows/Tests/badge.svg)

then these validation checks are working for the currently released version.

## PyCO2SYS calculations

All of the equilibrium constants and salt concentrations estimated from salinity are determined with [PyCO2SYS v1.8.1](https://PyCO2SYS.readthedocs.io/en/latest/) and therefore guaranteed to be identical to the values used there.  See the peer-reviewed [PyCO2SYS paper](https://doi.org/10.5194/gmd-15-15-2022) and its [validation docs](https://pyco2sys.readthedocs.io/en/latest/validate/) for more information on how these themselves are validated.
