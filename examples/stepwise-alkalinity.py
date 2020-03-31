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
concs = {
    'totalCarbonate': 1969.5e-6,
    'totalPhosphate': 0.1e-6,
    'totalSilicate': 0.7e-6,
    'totalAmmonia': 0,
    'totalH2Sulfide': 0,
}
volSample = 96
WhichKs = 10
WhoseKSO4 = 1
WhoseKF = 1
WhoseTB = 2

# Simulate junk first just to get correctly sized arrays
volAcid, emf, tempK  = calk.simulate.titration(acidVolStep=acidVolStep,
    alk0=alk0, buretteCorrection=buretteCorrection, concAcid=concAcid,
    emf0=emf0, maxVolAcid=maxVolAcid, pSal=pSal, tempK=tempK,
    volSample=volSample, **concs)

this = calk.Potentiometric(volAcid, emf, tempK)
this.add(pSal=pSal, volSample=volSample)
this.get_concTotals(**concs, WhichKs=WhichKs, WhoseTB=WhoseTB)
this.get_eqConstants(WhichKs=WhichKs, WhoseKSO4=WhoseKSO4, WhoseKF=WhoseKF)
this.solve_alkalinity(concAcid=concAcid)
this.get_alkSteps()
# todo: add other inputs (e.g. DIC)
# todo: do all checks for stepwiseAlkalinity

# Alkalinity anomaly stuff
massAcid = calk.convert.vol2massAcid(volAcid, tempK,
    buretteCorrection=buretteCorrection)
massSample = calk.convert.vol2massSample(volSample, tempK[0], pSal)
h = calk.solve.emf2h(emf, emf0, tempK)
mu = calk.solve.mu(massAcid, massSample)
concTotals = calk.concentrations.concTotals(pSal, **concs,
    WhichKs=WhichKs, WhoseTB=WhoseTB)
eqConstants = calk.dissociation.eqConstants(tempK, pSal, concTotals,
    WhichKs=WhichKs, WhoseKSO4=WhoseKSO4, WhoseKF=WhoseKF)
alkSim = calk.simulate.alk(h, mu, concTotals, eqConstants)[0]
alk0Sim = (alkSim + massAcid*concAcid/(massAcid + massSample))/mu

# Make the stepwise alkalinity plot
fig, ax = plt.subplots()
ax.scatter(volAcid, alk0Sim*1e6)
