# Test all solvers with a simulated titration

**Code:** [compare-all-solvers.py](https://github.com/mvdh7/calkulate/blob/master/examples/compare-all-solvers.py)

Here we import a titration simulated by [D81](../../references/#D81) and use it to try out all of Calkulate's alkalinity solvers.

Begin by importing the simulated titration dataset (this is hard-coded into Calkulate, so no separate files are required):

```python
import numpy as np
import calkulate as calk

# Import D81 simulated titration with phosphate
massAcid, pH, tempK, massSample, concAcid, pSal, alk, concTotals, eqConstants \
    = calk.io.Dickson1981(withPhosphate=False)
```

The dataset contains pH directly, rather than the measured EMF that Calkulate's solvers use. So let's convert the pH into EMF using a made-up EMF° value of 660 mV:

```python
# Define EMF0 & convert pH to EMF, as if we'd done a potentiometric titration
emf0 = 660.0
h = 10.0**-pH
emf = calk.solve.h2emf(h, emf0, tempK)
```

Now we can run all of the solvers on the simulated titration. As the solver functions all take the same set of inputs, we can put the inputs into the tuple `solveArgs` and splat them into the functions to save on writing them all out four times:

```python
# Solve for total alkalinity etc. with every solver and the true concAcid
solveArgs = (massAcid, emf, tempK, massSample, concAcid, concTotals,
    eqConstants)
alk_complete, emf0_complete = calk.solve.complete(*solveArgs)['x']
alk_DAA03, f_DAA03 = calk.solve.DAA03(*solveArgs)['x']
alk_Dickson1981, totalCarbonate_Dickson1981, f_Dickson1981 = \
    calk.solve.Dickson1981(*solveArgs)['x']
alk_halfGran, emf0_halfGran = calk.solve.halfGran(*solveArgs)['x']
```

Printing out the results, it becomes immediately apparent that the solvers do not all return the same alkalinity, and only the complete calculation method returns the *correct* alkalinity:

```python
# Print out results nicely
print('Total alkalinity in micromol/kg-sw:')
print(('{:^11} '*5).format('True', 'Complete', 'DAA03', 'Dickson1981',
    'halfGran'))
print('With true concAcid:')
print(('{:^11.2f} '*5).format(*np.array([alk, alk_complete, alk_DAA03,
    alk_Dickson1981, alk_halfGran])*1e6))
```

Why does this occur? The differences are mostly due to the different sets of equilibria that are accounted for in each solver. The Dickson CRM method [[DAA03](../../references/#DAA03)] appears to be particularly bad because it assumes that all CO<sub>2</sub> has been degassed from the sample, while the simulation is for a closed-cell titration with no CO<sub>2</sub> loss.

The first calculation above was done using the acid concentration that was declared by [D81](../../references/#D81) for the simulation. We could also compare the calibration solvers:

```python
# Solve for the acid concentration with every solver
calibrateArgs = (massAcid, emf, tempK, massSample, alk, concTotals,
    eqConstants)
concAcid_complete = calk.calibrate.concAcid(*calibrateArgs,
    solver='complete')['x'][0]
concAcid_DAA03 = calk.calibrate.concAcid(*calibrateArgs,
    solver='DAA03')['x'][0]
concAcid_Dickson1981 = calk.calibrate.concAcid(*calibrateArgs,
    solver='Dickson1981')['x'][0]
concAcid_halfGran = calk.calibrate.concAcid(*calibrateArgs,
    solver='halfGran')['x'][0]
```

Each solver returns a different `concAcid` value, for the same reason as why each solver returns a different alkalinity when a uniform `concAcid` input is used. However, if we now solve for alkalinity again with each solver but using the appropriate `concAcid` to each one:

```python
# Solve for total alkalinity etc. with every solver and its own concAcid
alk_complete_cal, emf0_complete_cal = calk.solve.complete(massAcid, emf, tempK,
    massSample, concAcid_complete, concTotals, eqConstants)['x']
alk_DAA03_cal, f_DAA03_cal = calk.solve.DAA03(massAcid, emf, tempK, massSample,
    concAcid_DAA03, concTotals, eqConstants)['x']
alk_Dickson1981_cal, totalCarbonate_Dickson1981_cal, f_Dickson1981_cal = \
    calk.solve.Dickson1981(massAcid, emf, tempK, massSample,
    concAcid_Dickson1981, concTotals, eqConstants)['x']
alk_halfGran_cal, emf0_halfGran_cal = calk.solve.halfGran(massAcid, emf, tempK,
    massSample, concAcid_halfGran, concTotals, eqConstants)['x']
```

We now find that the self-calibrated solvers all return the correct alkalinity (i.e. the value that the simulation was based upon):

```python
# Print out self-calibrated results
print('With self-calibrated concAcid:')
print(('{:^11.2f} '*5).format(*np.array([alk, alk_complete_cal, alk_DAA03_cal,
    alk_Dickson1981_cal, alk_halfGran_cal])*1e6))
```

The moral of the story, which is echoed in some other aspects of titration data processing such as volume calibrations, is that using a consistent approach for calibration and sample processing is more important for getting a good alkalinity value than exactly which calculation approach is used.
