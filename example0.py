import numpy as np
import calkulate as calk

# Import D81 simulated titration with phosphate
massAcid, pH, tempK, massSample, concAcid, pSal, alk, concTotals, eqConstants \
    = calk.io.Dickson1981(withPhosphate=False)
totalCarbonate = concTotals['C']

# Define EMF0 & convert pH to EMF, as if we'd done a potentiometric titration
emf0 = 660.0
h = 10.0**-pH
emf = calk.solve.h2emf(h, emf0, tempK)

# Solve for total alkalinity etc. with every solver and the true concAcid
solveArgs = (massAcid, emf, tempK, massSample, concAcid, concTotals,
    eqConstants)
alk_complete, emf0_complete = calk.solve.complete(*solveArgs)['x']
alk_DAA03, f_DAA03 = calk.solve.DAA03(*solveArgs)['x']
alk_Dickson1981, totalCarbonate_Dickson1981, f_Dickson1981 = \
    calk.solve.Dickson1981(*solveArgs)['x']
alk_halfGran, emf0_halfGran = calk.solve.halfGran(*solveArgs)['x']

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

# Print out results nicely
print('Total alkalinity in micromol/kg-sw:')
print(('{:^11} '*5).format('True', 'Complete', 'DAA03', 'Dickson1981',
    'halfGran'))
print('With true concAcid:')
print(('{:^11.2f} '*5).format(*np.array([alk, alk_complete, alk_DAA03,
    alk_Dickson1981, alk_halfGran])*1e6))
print('With self-calibrated concAcid:')
print(('{:^11.2f} '*5).format(*np.array([alk, alk_complete_cal, alk_DAA03_cal,
    alk_Dickson1981_cal, alk_halfGran_cal])*1e6))
