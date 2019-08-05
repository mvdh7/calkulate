import numpy as np
from matplotlib import pyplot as plt
import calkulate as calk

# Import D81 simulated titration with and without phosphate
(massAcid0, pH0, tempK0, massSample0, concAcid0, pSal0, alk0, concTotals0,
    eqConstants0) = calk.io.Dickson1981(withPhosphate=False)
(massAcid1, pH1, tempK1, massSample1, concAcid1, pSal1, alk1, concTotals1,
    eqConstants1) = calk.io.Dickson1981(withPhosphate=True)

# Simulate it again with Calkulate and round to the same precision as D81
pHSim0 = calk.simulate.pH(massAcid0, massSample0, concAcid0, alk0, concTotals0,
    eqConstants0)
pHSim1 = calk.simulate.pH(massAcid1, massSample1, concAcid1, alk1, concTotals1,
    eqConstants1)
pHSim0 = np.round(pHSim0, decimals=6)
pHSim1 = np.round(pHSim1, decimals=6)

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
