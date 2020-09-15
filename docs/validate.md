# Validating Calkulate

## Automatic tests

Several of the [automatic tests](https://github.com/mvdh7/calkulate/tree/master/tests) include checks on whether Calkulate does self-consistent calibrations, for example solving a sample where the titrant was calibrated with itself, which should return exactly the correct alkalinity.

The tests also include a check that Calkulate can return the correct alkalinity and titrant molinity from the two simulated titration datasets published by [D81](../references/#d).

So if the build badge below says "passing":

[![Build Status](https://travis-ci.org/mvdh7/calkulate.svg?branch=master)](https://travis-ci.org/mvdh7/calkulate)

then these validation checks are working for the currently released version.

## PyCO2SYS calculations

All of the equilibrium constants and salt concentrations estimated from salinity are determined with [PyCO2SYS](https://PyCO2SYS.rtfd.io) and therefore guaranteed to be identical to the values used there.  See the [PyCO2SYS validation docs](https://pyco2sys.readthedocs.io/en/latest/validate/) for more information on how these themselves are validated.
