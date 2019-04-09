import calkulate as calk

# Import Dickson (1981) dataset
Macid, pH, tempK, Msamp, Cacid, psal, AT_cert, XT, KXF = calk.io.Dickson1981()

# Convert pH to EMF, as if it were a potentiometric titration
h = 10.0**-pH
emf0 = 660.0
emf = calk.solve.h2emf(h, emf0, tempK)

# Recover the acid concentration, as if calibrating
Cacid_cal = calk.calibrate.complete(Macid, emf, tempK, Msamp, AT_cert,
    XT, KXF)['x'][0]

# Recover the total alkalinity and EMF0 from the titration data
AT_cal, emf0_cal = calk.solve.complete(Macid, emf, tempK, Msamp, Cacid, 
    XT, KXF)['x']

# Re-simulate the Dickson (1981) dataset
pH_sim = calk.simulate.pH(Macid, Msamp, Cacid, AT_cert, XT, KXF)
