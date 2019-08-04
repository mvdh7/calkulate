# Calibrating the acid titrant

The acid titrant can be calibrated using measurements of a sample with known alkalinity: we simply solve for the acid concentration that produces the known alkalinity value.

In Calkulate this can be done using the acid concentration calibration function:

```python
concAcid = calk.calibrate.concAcid(massAcid, emf, tempK, massSample, alkCert,
    concTotals, eqConstants, solver='complete')['x'][0]
```

As with the [solver functions](../solvers), the entire `OptimizeResult` from [SciPy's least-squares optimisation function](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html) is returned, so the actual acid concentration value is stored in the field `'x'`.

By default this calibration function uses the [complete calculation](../solvers/#complete-complete-calculation) method, but you can use any of the other solvers by changing the `solver` input. The (case-insensitive) options are:

  * `solver='complete'`: [complete calculation](../solvers/#complete-complete-calculation) (default);
  * `solver='DAA03'`: [Dickson CRM method](../solvers/#daa03-dickson-crm-method);
  * `solver='Dickson1981'`: [closed-cell titrations](../solvers/#dickson1981-closed-cell-titrations);
  * `solver='halfGran'`: [half-Gran method](../solvers/#halfgran-half-gran-method).
