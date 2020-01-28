import numpy as np
import matplotlib.pyplot as plt
import calkulate as calk

# Simulate a VINDTA-style titration
acidVolStep = 0.15
alk0 = 2238.6e-6
buretteCorrection = 1
concAcid = 0.1
emf0 = 660
maxVolAcid = 4.1
pSal = 33.571
tempK = 298.15
totalCarbonate = 2031.53e-6
totalPhosphate = 0.31e-6
totalSilicate = 2.5e-6
volSample = 100
volAcid, emf, tempK  = calk.simulate.titration(acidVolStep=acidVolStep,
    alk0=alk0, buretteCorrection=buretteCorrection, concAcid=concAcid,
    emf0=emf0, maxVolAcid=maxVolAcid, pSal=pSal, tempK=tempK,
    totalCarbonate=totalCarbonate, totalPhosphate=totalPhosphate,
    totalSilicate=totalSilicate, volSample=volSample)

# Now replace results with totalCarbonate as array
totalCarbonateVec = np.full_like(emf, totalCarbonate) # flat
totalCarbonateVec = totalCarbonate*(0.98 + np.exp(-volAcid)*0.02)
volAcid, emf, tempK = calk.simulate.titration(acidVolStep=acidVolStep,
    alk0=alk0, buretteCorrection=buretteCorrection, concAcid=concAcid,
    emf0=emf0, maxVolAcid=maxVolAcid, pSal=pSal, tempK=tempK,
    totalCarbonate=totalCarbonateVec, totalPhosphate=totalPhosphate,
    totalSilicate=totalSilicate, volSample=volSample)

# Solve it for alkalinity
massSample = volSample*calk.density.sw(tempK[0], pSal)*1e-3
massAcid = buretteCorrection*volAcid*calk.density.acid(tempK)*1e-3
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
fig, ax = plt.subplots(2, 1)
ax[0].scatter(massAcid, alk0Est*1e6)
ax[0].set_xlim([0, np.max(massAcid)])
ax[0].plot([0, np.max(massAcid)], alk0Solved*np.array([1, 1])*1e6, c='k')
ax[1].scatter(massAcid, totalCarbonateVec*1e6)
ax[1].set_xlim([0, np.max(massAcid)])
ax[1].plot([0, np.max(massAcid)], totalCarbonate*np.array([1, 1])*1e6, c='k')
