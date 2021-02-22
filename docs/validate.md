# Validating Calkulate

## Automatic tests

Several of the [automatic tests](https://github.com/mvdh7/calkulate/tree/master/tests) include checks on whether Calkulate does self-consistent calibrations, for example solving a sample where the titrant was calibrated with itself, which should return exactly the correct alkalinity.

So if the Tests badge below says "passing":

![Tests](https://github.com/mvdh7/calkulate/workflows/Tests/badge.svg)

then these validation checks are working for the currently released version.

## PyCO2SYS calculations

All of the equilibrium constants and salt concentrations estimated from salinity are determined with [PyCO2SYS v1.6.0](https://PyCO2SYS.readthedocs.io/en/latest/) and therefore guaranteed to be identical to the values used there.  See the [PyCO2SYS validation docs](https://pyco2sys.readthedocs.io/en/latest/validate/) for more information on how these themselves are validated.
