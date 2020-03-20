import calkulate as calk

# Certified reference material (CRM) batch 144 values from:
#     https://www.nodc.noaa.gov/ocads/oceans/Dickson_CRM/batches.html
pSal = 33.571
totalCarbonate = 2031.53e-6
totalPhosphate = 0.31e-6
totalSilicate = 2.5e-6
alkCert = 2238.60e-6
totalAmmonia = 0 # 1.0e-6 # use non-zero value just for testing
totalH2Sulfide = 0 # 1.0e-6 # use non-zero value just for testing

# Import a VINDTA-generated .dat file from a real CRM-144 analysis
volSample = 99.981 # ml
datFile = 'datfiles/CRM-144-0435-4.dat'
massAcid, emf, tempK, massSample, concTotals, eqConstants = \
    calk.vindta.prep(datFile, volSample, pSal, totalCarbonate, totalPhosphate,
    totalSilicate, totalAmmonia=totalAmmonia, totalH2Sulfide=totalH2Sulfide)

# Calibrate the acid concentration based on this CRM measurement
concAcid = calk.calibrate.concAcid(massAcid, emf, tempK, massSample, alkCert,
    concTotals, eqConstants)['x'][0]

# Check that the calibrated acid concentration returns the correct alkalinity
alk = calk.solve.complete(massAcid, emf, tempK, massSample, concAcid,
    concTotals, eqConstants)['x'][0]
print('Certified alkalinity  = {:.2f} micromol/kg-sw'.format(alkCert*1e6))
print('Calibrated alkalinity = {:.2f} micromol/kg-sw'.format(alk*1e6))

# Visualise "everything" about this titration
calk.plot.everything(datFile, volSample, pSal, totalCarbonate, totalPhosphate,
    totalSilicate, concAcid, totalAmmonia=totalAmmonia,
    totalH2Sulfide=totalH2Sulfide)
