import calkulate as calk
import numpy as np

# Do a real VINDTA .dat file
datfile = 'datfiles/0-0  0  (0)CRM-144-0435-4.dat'
#datfile = 'datfiles/calk_simtit.dat'
#datfile = 'Ruth Pre-Locate/2019_0305_67_0_0_1.dat'
Vacid, emf, Tk = calk.io.vindta(datfile)
Dacid = calk.density.acid(Tk[0]) # kg/l
Macid = Vacid * Dacid * 1e-3 # kg
Vsamp = 100. # ml
psal = 33.571
Msamp = Vsamp * calk.density.sw(Tk[0], psal) / 1e6 # kg

# Get *XT and *KX
AT_cert = 0.00223860
# ^^^ solver works with other values but not with the real one???
CT  = 0.00203153
PT  = 3.1e-7
SiT = 2.5e-6
XT = calk.concentrations.XT(psal, CT, PT, SiT)
KXF = calk.dissociation.KXF(Tk, psal, XT)

Macid, emf, tempK, Msamp, XT, KXF = calk.vindta.prep(datfile, Vsamp, psal, 
    CT, PT, SiT)

Cacid = 0.1

ATg, emf0g, _, pHg = calk.solve.guessGran(Macid, emf, tempK, Msamp, Cacid)
L = np.logical_and(pHg > 3, pHg < 4)
KXL = {k: v[L] for k, v in KXF.items()}

Cacid = calk.calibrate.halfGran(Macid, emf, tempK, Msamp, AT_cert,
    XT, KXF)['x'][0]
test1 = calk.solve.halfGran(Macid, emf, tempK, Msamp, Cacid, XT, KXF)


## Solve
#
#Cacid = calk.vindta.completeCRM(datfile, Vsamp, AT_cert, psal, CT, PT, SiT)[0]
#test = calk.vindta.complete(datfile, Vsamp, Cacid, psal, CT, PT, SiT)['x']
#H = calk.solve.emf2h(emf, test[1], Tk)
#simAT = calk.vindta.simAT(Macid, Tk, H, Msamp, psal, CT, PT, SiT)[0]
