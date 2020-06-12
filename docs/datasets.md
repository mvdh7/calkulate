# Titration datasets

## Import with Calkulate

If your titration dataset contains all the necessary information described in [Input data formatting](../inputs) in correctly labelled columns, then using Calkulate is straightforward:

    :::python
    import calkulate as calk

    # If your titration table is in a file, just provide the path and name:
    tdata = calk.Dataset("path/to/titration_table.csv")

    # If you have already imported your titration table as a pandas DataFrame
    # or anything that can be converted into one, just provide the table itself:
    tdata = calk.Dataset(titration_table)

    # The only Calkulate command you may ever need:
    tdata.calibrate_and_solve()

The result `tdata` consists of three components:

  1.  The titration table is now stored as a standard [pandas DataFrame](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html) in `tdata.table`.
  2.  All of the individual titration data files have also been imported and are stored in `tdata.titrations` (see [Individual titrations](../titrations) for more details).
  3.  A second DataFrame containing metadata about each analysis batch is stored in `tdata.batches`.

Read on for more detail on the calibration and alkalinity-solving step.

## Calibrate the titrant molinities

Assuming that the titration table contains some rows with a certified alkalinity value (reference material) and others without (samples), we can now

  1.  Calibrate the titrant molinity for each indiviual reference,
  2.  Calculate the average titrant molinity across each analysis batch, and
  3.  Add the appropriate batch-average titrant molinities into a `titrant_molinity` column in the titration table

with the command:

    :::python
    tdata.calibrate()

## Solve for alkalinity

To then solve every titration in the table for its total alkalinity using these calibrated titrant molinities and store the result in a `alkalinity` column in the titration table:

    :::python
    tdata.solve()
    