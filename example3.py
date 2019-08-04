import pandas as pd
import numpy as np
import calkulate as calk

# Import dataset summary, and set titration subsample volume and solver method
dataset = pd.read_csv('datfiles/exampleDataset.csv')
volSample = 99.956 # ml
method = 'complete'

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

# Solve every titration for alkalinity
dataset['alk'] = np.nan
for i in dataset.index:
    dataset.loc[i, 'alk'] = calk.vindta.alk(
        'datfiles/{}.dat'.format(dataset.fileName[i]), volSample,
        dataset.concAcidMean[i], dataset.pSal[i],
        dataset.totalCarbonate[i]*1e-6, dataset.totalPhosphate[i]*1e-6,
        dataset.totalSilicate[i]*1e-6, solver=method)['x'][0]
dataset.loc[:, 'alk'] *= 1e6 # convert to micromol/kg
