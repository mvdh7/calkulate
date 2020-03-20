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

---

## Wrapper for .dat file inputs

If your data is formatted as a [VINDTA-style .dat file](../io/#what-format-is-that), you can more easily implement any of the solvers above using the wrapped `concAcid` function in the `datfile` module:

```python
concAcidOptResult = calk.datfile.alk(datFile, volSample, alkCert, pSal, totalCarbonate,
    totalPhosphate, totalSilicate, solver='complete', buretteCorrection=1, tempKForce=None)
```

Here, the input `datFile` is the file name (and path to it, if necessary) of the titration data file that you want to solve. Just like in the [calibration function](#calibrating-the-acid-titrant), the case-insensitive optional `solver` input (defaults to the [complete calculation](../solvers/#complete-complete-calculation) method) can be set as the name of any of the methods available in Calkulate.

The output `concAcidOptResult` contains the full optimisation output that would be obtained by using the `calibrate.concAcid` method described above.
