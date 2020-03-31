import numpy as np
import matplotlib.pyplot as plt
import calkulate as calk

# Simulate a VINDTA-style titration
acidVolStep = 0.15
alk0 = 2279.7e-6
buretteCorrection = 1
concAcid = 0.1
emf0 = 638.5
maxVolAcid = 4.1
pSal = 33.571
tempK = 298.15
totalCarbonate = 1969.5e-6
totalPhosphate = 0.1e-6
totalSilicate = 0.7e-6
volSample = 96

# Simulate junk first just to get correctly sized arrays
volAcid0, emf, tempK  = calk.simulate.titration(acidVolStep=acidVolStep,
    alk0=alk0, buretteCorrection=buretteCorrection, concAcid=concAcid,
    emf0=emf0, maxVolAcid=maxVolAcid, pSal=pSal, tempK=tempK,
    totalCarbonate=totalCarbonate, totalPhosphate=totalPhosphate,
    totalSilicate=totalSilicate, volSample=volSample)

# Now replace results with totalCarbonate as array and simulate again
totalCarbonateVec = np.full_like(emf, totalCarbonate) # flat
# totalCarbonateVec = totalCarbonate*(volSample - volAcid0)/volSample
# totalCarbonateVec = totalCarbonate*(0.98 + np.exp(-volAcid0)*0.02)
# totalCarbonateVec = totalCarbonate*(0.4 + np.exp(-volAcid0*2)*0.3)
acidMult = 1.05
volAcid, emf, tempK = calk.simulate.titration(acidVolStep=acidVolStep*acidMult,
    alk0=alk0, buretteCorrection=buretteCorrection, concAcid=concAcid,
    emf0=emf0, maxVolAcid=maxVolAcid*acidMult, pSal=pSal, tempK=tempK,
    totalCarbonate=totalCarbonateVec, totalPhosphate=totalPhosphate,
    totalSilicate=totalSilicate*0, volSample=volSample, extraVolAcid=0.0)

# Solve it for alkalinity
massSample = volSample*calk.density.sw(tempK[0], pSal)*1e-3
massAcid = buretteCorrection*volAcid0*calk.density.acid(tempK)*1e-3
concTotals = calk.concentrations.concTotals(pSal, totalCarbonate,
    totalPhosphate, totalSilicate)
eqConstants = calk.dissociation.eqConstants(tempK, pSal, concTotals)
alk0Solved, emf0Solved = calk.solve.complete(massAcid, emf, tempK, massSample,
    concAcid, concTotals, eqConstants)['x']

# Estimate alkalinity at every titration point
hSolved = calk.solve.emf2h(emf, emf0Solved, tempK)
pHSolved = -np.log10(hSolved)
mu = calk.solve.mu(massAcid, massSample)
alkSim, alkSimComponents = calk.simulate.alk(hSolved, mu, concTotals,
    eqConstants)
alk0Est = alkSim/mu + massAcid*concAcid/massSample

# Plot differences
print('True alkalinity = {:.2f}'.format(alk0*1e6))
print('Calk alkalinity = {:.2f}'.format(alk0Solved*1e6))
fig, ax = plt.subplots(2, 1, figsize=(8, 13))
ax[0].scatter(massAcid*1e3, alk0Est*1e6)
ax[0].set_xlim([0, np.max(massAcid)*1e3])
ax[0].plot([0, np.max(massAcid)*1e3], alk0Solved*np.array([1, 1])*1e6, c='r')
ax[0].set_xlabel('Mass acid / g')
ax[1].scatter(massAcid*1e3, totalCarbonateVec*1e6)
ax[1].set_xlim([0, np.max(massAcid*1e3)])
ax[1].plot([0, np.max(massAcid*1e3)], totalCarbonate*np.array([1, 1])*1e6,
    c='r')
ax[1].set_xlabel('Mass acid / g')
