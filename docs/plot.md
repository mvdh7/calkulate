# Visualising the results

After running `calkulate()` (or `calibrate()` and `solve()`) on your data, Calkulate contains some plotting functions to help visualise the results.  More will be added in time, and contributions are welcome!

## Calibrated `titrant_molinity`

```python
fig, ax = calk.plot.titrant_molinity(
    data, xvar=None, show_bad=True, show_batches=True, figure_fname=None,
)
```

The required argument `data` is the [metadata table](../metadata) as a pandas DataFrame or Calkulate Dataset.

Optional inputs:

  * `xvar`: name of column to use as the x-axis variable.
  * `show_bad`: whether or not to show values where `data.reference_good == False`.
  * `show_batches`: whether or not to show batch-averaged `titrant_molinity` values.
  * `figure_fname`: if provided, save figure to this filename.

## Measured − certified `alkalinity_offset`

```python
fig, ax = calk.plot.alkalinity_offset(
    data, xvar=None, show_bad=True, show_batches=True, figure_fname=None,
)
```

The required argument `data` is the [metadata table](../metadata) as a pandas DataFrame or Calkulate Dataset.

Optional inputs:

  * `xvar`: name of column to use as the x-axis variable.
  * `show_bad`: whether or not to show values where `data.reference_good == False`.
  * `show_batches`: whether or not to show batch-averaged `titrant_molinity` values.
  * `figure_fname`: if provided, save figure to this filename.


