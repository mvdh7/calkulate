# Simulate a titration with Calkulate

**Code:** [resimulate-d81.py](https://github.com/mvdh7/calkulate/blob/master/examples/resimulate-d81.py)

Here we use Calkulate to re-simulate the titration data reported by [D81](../../references/#D81).

First, we import the [D81](../../references/#D81) simulation datasets with and without phosphate (these are hard-coded into Calkulate, so no separate files are required):

```python
import numpy as np
from matplotlib import pyplot as plt
import calkulate as calk

# Import D81 simulated titration with and without phosphate
(massAcid0, pH0, tempK0, massSample0, concAcid0, pSal0, alk0, concTotals0,
    eqConstants0) = calk.io.Dickson1981(withPhosphate=False)
(massAcid1, pH1, tempK1, massSample1, concAcid1, pSal1, alk1, concTotals1,
    eqConstants1) = calk.io.Dickson1981(withPhosphate=True)
```

Now, we can simulate the datasets again with Calkulate. We round the results to the same number of significant figures as [D81](../../references/#D81), to make it easier to check how consistent the results are:

```python
# Simulate it again with Calkulate and round to the same precision as D81
pHSim0 = calk.simulate.pH(massAcid0, massSample0, concAcid0, alk0, concTotals0,
    eqConstants0)
pHSim1 = calk.simulate.pH(massAcid1, massSample1, concAcid1, alk1, concTotals1,
    eqConstants1)
pHSim0 = np.round(pHSim0, decimals=6)
pHSim1 = np.round(pHSim1, decimals=6)
```

We visualise the difference between our new simulation and the original pH data:

```python
# Plot the differences between Calkulate and D81 simulations
fig, ax = plt.subplots()
ax.scatter(massAcid0*1e3, (pHSim0 - pH0)*1e3, marker='+',
    label='No phosphate')
ax.scatter(massAcid1*1e3, (pHSim1 - pH1)*1e3, marker='x',
    label='With phosphate')
ax.plot([0, np.max(massAcid0)*1e3], [0, 0], c='k')
ax.set_xlabel('Acid mass / g')
ax.set_ylabel('[pH(Calkulate) $-$ pH(D81)] × 10$^3$')
ax.legend()
```

The result:

<img src='../../img/resimulate-d81.png' />

All of the no-phosphate data simulated by Calkulate agree exactly with [D81](../../references/#D81).

Most of the with-phosphate data agree, but a couple of the points differ slightly. It's hard to imagine any explanation for these isolated discrepancies other than that there may be a couple of minor typographical errors in [D81](../../references/#D81)'s Table 4.
