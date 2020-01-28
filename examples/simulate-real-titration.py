from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
import calkulate as calk

# Read titration data from .dat file.
# It's a CRM: self-calibrate it following examples/self-calibrate-crm.py
# Certified reference material (CRM) batch 144 values from:
#     https://www.nodc.noaa.gov/ocads/oceans/Dickson_CRM/batches.html
pSal = 33.571
totalCarbonate = 2031.53e-6
totalPhosphate = 0.31e-6
totalSilicate = 2.5e-6
alkCert = 2238.60e-6
volSample = 99.981 # ml
datFile = 'datfiles/CRM-144-0435-4.dat'
massAcid, emf, tempK, massSample, concTotals, eqConstants = \
    calk.vindta.prep(datFile, volSample, pSal, totalCarbonate, totalPhosphate,
    totalSilicate)
concAcid = calk.calibrate.concAcid(massAcid, emf, tempK, massSample, alkCert,
    concTotals, eqConstants)['x'][0]
alkCheck, emf0Check = calk.solve.complete(massAcid, emf, tempK, massSample,
    concAcid, concTotals, eqConstants)['x']

test = calk.simulate.titration()
print(test)

# #%% Now, simulate this titration
# mu = calk.solve.mu(massAcid, massSample)
# # Set simulation conditions
# concTotalsSimulator = deepcopy(concTotals)
# # concTotalsSimulator['Si'] = 10e-6
# concTotalsSimulator['C'] = np.full_like(eqConstants['C1'], 2031.53e-6)*mu**12
# # Set solving conditions
# concTotalsSolver = deepcopy(concTotals)
# # Simulate and solve
# pHSim = calk.simulate.pH(massAcid, massSample, concAcid, alkCert,
#     concTotalsSimulator, eqConstants)
# hSim = 10.0**-pHSim
# emfSim = calk.solve.h2emf(hSim, emf0, tempK)
# alk0SimSolved, emf0SimSolved = calk.solve.complete(massAcid, emfSim, tempK,
#     massSample, concAcid, concTotalsSolver, eqConstants)['x']
# hSimSolved = calk.solve.emf2h(emfSim, emf0SimSolved, tempK)
# alkTitSimSolved = calk.simulate.alk(hSimSolved, mu, concTotalsSolver,
#     eqConstants)[0]
# alkEst = (alkTitSimSolved + massAcid*concAcid/(massAcid + massSample))/mu
# # Plot difference
# fig, ax = plt.subplots(2, 1)
# ax[0].scatter(massAcid, alkEst*1e6)
# ax[0].set_xlim([0, np.max(massAcid)])
# ax[0].plot([0, np.max(massAcid)], alk0SimSolved*np.array([1, 1])*1e6, c='k')

