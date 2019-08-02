# Solving for alkalinity

Several different approaches for determining total alkalinity from titration data are implemented in Calkulate. They are implemented by functions in Calkulate's `solve` module, and they all share a common syntax:

```python
OptimizeResult = calk.solve.method(massAcid, emf, tempK, massSample, concAcid, XT, KXF)
```

The formats of the input variables are described in the documentation section on [variables and conventions](../conventions).

The output `OptimizeResult` is usually (except for the `halfGran` method, see later) the output of SciPy's `optimize.least_squares` function. As described in [the SciPy documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html), the actual fitted alkalinity value (and any auxiliary fitted parameters, which vary by method) is contained within the field `OptimizeResult['x']`. The alkalinity is always the first item in this list, i.e. `OptimizeResult['x'][0]`.
