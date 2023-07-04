# Visualising the results

After running `calkulate()` (or `calibrate()` and `solve()`) on your data, Calkulate contains some plotting functions to help visualise the results.  More will be added in time, and contributions are welcome!

## Dataset plots

### Calibrated `titrant_molinity`

```python
fig, ax = calk.plot.titrant_molinity(
    ds, xvar=None, show_bad=True, show_batches=True, figure_fname=None,
)
```

The required argument `ds` is the [metadata table](../metadata) as a pandas DataFrame or Calkulate Dataset.

Optional inputs:

  * `xvar`: name of column to use as the x-axis variable.
  * `show_bad`: whether or not to show values where `ds.reference_good == False`.
  * `show_batches`: whether or not to show batch-averaged `titrant_molinity` values.
  * `figure_fname`: if provided, save figure to this filename.

### Measured − certified `alkalinity_offset`

```python
fig, ax = calk.plot.alkalinity_offset(
    ds, xvar=None, show_bad=True, show_batches=True, figure_fname=None,
)
```

The required argument `ds` is the [metadata table](../metadata) as a pandas DataFrame or Calkulate Dataset.

Optional inputs:

  * `xvar`: name of column to use as the x-axis variable.
  * `show_bad`: whether or not to show values where `ds.reference_good == False`.
  * `show_batches`: whether or not to show batch-averaged `titrant_molinity` values.
  * `figure_fname`: if provided, save figure to this filename.

## Titration plots

To investigate an individual titration in more detail, first generate a `Titration` from the relevant row of your `Dataset`:

```python
tt = ds.to_Titration(index)
```

where `index` is the index value for the row you are interested in.

A series of figures can then be plotted for the titration in question:

  * `tt.plot_emf()`: how EMF changes through the titration.
  * `tt.plot_pH()`: how pH changes through the titration.
  * `tt.plot_gran_alkalinity()`: the Gran-plot initial alkalinity estimate.
  * `tt.plot_gran_emf0()`: the Gran-plot initial EMF<sup>0</sup> estimate.
  * `tt.plot_alkalinity()`: the total alkalinity calculated from each titration data point.
  * `tt.plot_components()`: how every equilibrating component of the solution changes throughout the titration.
