import matplotlib.pyplot as plt
import calkulate as calk

# Define inputs to simulate a potentiometric titration
acidVolStep = 0.15
alk0 = 2279.7e-6
buretteCorrection = 1
concAcid = 0.1
emf0 = 638.5
maxVolAcid = 4.1
pSal = 33.571
tempK = 298.15
volSample = 96
concs = {
    'totalCarbonate': 1969.5e-6,
    'totalPhosphate': 0.1e-6,
    'totalSilicate': 0.7e-6,
    'totalAmmonia': 0.1e-6,
    'totalH2Sulfide': 0.1e-6,
}
Ks = {'WhichKs': 10, 'WhoseKSO4': 1, 'WhoseKF': 1, 'WhoseTB': 2}

# Simulate the arrays then put them into a Potentiometric object
volAcid, emf, tempK = calk.simulate.titration(acidVolStep=acidVolStep,
    alk0=alk0, buretteCorrection=buretteCorrection, concAcid=concAcid,
    emf0=emf0, maxVolAcid=maxVolAcid, pSal=pSal, tempK=tempK,
    volSample=volSample, **concs, **Ks)
this = calk.titration.Potentiometric(volAcid, emf, tempK, pSal, volSample,
                                      **concs, **Ks)

# # Alternatively, just simulate directly into a Potentiometric object
# this = calk.simulate.Potentiometric(acidVolStep=acidVolStep,
#     alk0=alk0, buretteCorrection=buretteCorrection, concAcid=concAcid,
#     emf0=emf0, maxVolAcid=maxVolAcid, pSal=pSal, tempK=tempK,
#     volSample=volSample, **concs, **Ks)

solver = 'complete' # this is the only valid option for now
this.solve_alk(concAcid, solver=solver)
this.get_alkSteps(solver=solver)

# Plot the components of alkalinity along the way
fig, ax = plt.subplots(figsize=(5, 8))
calk.plot.alkComponents(this)

# Import titration from .dat file into a Potentiometric object
that = calk.datfile.Potentiometric('datfiles/CRM-144-0435-4.dat',
                                   volSample, pSal, **concs, **Ks)
that.calibrate_concAcid(2250.0e-6)
that.solve_alk(that.concAcidCert['complete'])
that.get_alkSteps()

# Plot the stepwise alkalinity estimates
fig, ax = plt.subplots(figsize=(5, 4))
calk.plot.alkSteps(that)
