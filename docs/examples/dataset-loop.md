# Calibrate and solve a series of measurements

**Code:** [dataset-loop.py](https://github.com/mvdh7/calkulate/blob/master/examples/dataset-loop.py) / **Data files:** [available from GitHub](https://github.com/mvdh7/calkulate/blob/master/datfiles)

Here we calibrate the acid concentration using data from several certified reference material (CRM) titrations, and then apply this to solve alkalinity for a series of real seawater samples.

It's helpful to create a summary spreadsheet that contains the relevant information about every sample and CRM, so we know which files to load and what to do with them. First, we import this summary, and set some general properties for the solver:

```python
import pandas as pd
import numpy as np
import calkulate as calk

# Import dataset summary, and set titration subsample volume and solver method
dataset = pd.read_csv('datfiles/exampleDataset.csv')
volSample = 99.956 # ml
method = 'complete'
```

Next, we loop through the CRMs, and calculate the best-fit acid concentration for each one. We could calculate the average value of these acid concentrations as our best estimate of its true value.

```python
# Get best-fit acid concentration for each CRM and then calculate average
dataset['concAcidPerCRM'] = np.nan
for i in dataset.index:
    if dataset.type[i] == 'CRM':
        dataset.loc[i, 'concAcidPerCRM'] = calk.vindta.concAcid(
            'datfiles/{}.dat'.format(dataset.fileName[i]),
            volSample, dataset.alkCert[i]*1e-6, dataset.pSal[i],
            dataset.totalCarbonate[i]*1e-6, dataset.totalPhosphate[i]*1e-6,
            dataset.totalSilicate[i]*1e-6, solver=method)['x'][0]
dataset['concAcidMean'] = dataset.concAcidPerCRM.mean()
```

Now that we have a calibrated acid concentration, we can loop through all of the titration data files and solve for alkalinity:

```python
# Solve every titration for alkalinity
dataset['alk'] = np.nan
for i in dataset.index:
    dataset.loc[i, 'alk'] = calk.vindta.alk(
        'datfiles/{}.dat'.format(dataset.fileName[i]), volSample,
        dataset.concAcidMean[i], dataset.pSal[i],
        dataset.totalCarbonate[i]*1e-6, dataset.totalPhosphate[i]*1e-6,
        dataset.totalSilicate[i]*1e-6, solver=method)['x'][0]
dataset.loc[:, 'alk'] *= 1e6 # convert to micromol/kg
```

Note that the CRMs have been calibrated using the average acid concentration, not using their own indiviual calibrations, so their alkalinities are not *exactly* the same as the certified values.

Finally, we can export the results back into a file, if required:

```python
# Put new results back into a CSV file
dataset.to_csv('datfiles/exampleDataset_calibrated.csv')
```
