import calkulate as calk

datFile = 'datfiles/CRM-144-0435-4.dat'
volSample = 100. # ml
pSal = 33.571
totalCarbonate = 0.00203153
totalPhosphate = 3.1e-7
totalSilicate = 2.5e-6
massAcid, emf, tempK, massSample, concTotals, eqConstants = \
    calk.vindta.prep(datFile, volSample, pSal, totalCarbonate, totalPhosphate,
    totalSilicate)
alkCert = 2230e-6
concAcid = calk.calibrate.concAcid(massAcid, emf, tempK, massSample, alkCert,
    concTotals, eqConstants, solver='halfGran')['x'][0]
calk.plot.everything(datFile, volSample, pSal, totalCarbonate, totalPhosphate,
    totalSilicate, concAcid)

#calk.io.writeDat('datfiles/test.dat', massAcid, emf, tempK)
