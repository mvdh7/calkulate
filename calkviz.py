import calkulate as calk
import numpy as np
from matplotlib import pyplot as plt

datfile = 'no-release/0-0  0  (0)CRM-144-0435-4.dat'
Vacid, emf, Tk = calk.io.vindta(datfile)
Dacid = calk.density.acid(Tk[0]) # kg/l
Macid = Vacid * Dacid * 1e-3 # kg
Vsamp = 100. # ml
psal = 33.571
Msamp = Vsamp * calk.density.sw(Tk[0], psal) / 1e6 # kg
AT_cert = 0.00223860
CT  = 0.00203153
PT  = 3.1e-7
SiT = 2.5e-6
XT = calk.concentrations.XT(psal, CT, PT, SiT)
KXF = calk.dissociation.KXF(Tk, psal, XT)
Macid, emf, tempK, Msamp, XT, KXF = calk.vindta.prep(datfile, Vsamp, psal, 
    CT, PT, SiT)
Cacid = calk.calibrate.complete(Macid, emf, tempK, Msamp, AT_cert, XT,
    KXF)['x'][0]
ATg, emf0g, _, pHg = calk.solve.guessGran(Macid, emf, tempK, Msamp, Cacid)
L = np.logical_and(pHg > 3, pHg < 4)
KXL = {k: v[L] for k, v in KXF.items()}
ATsolved = calk.solve.complete(Macid, emf, tempK, Msamp, Cacid, XT,
    KXF)['x'][0]
f1 = calk.solve.f1(Macid, emf, tempK, Msamp)

fig, ax = plt.subplots(3, 2)

ax[0, 0].scatter(Macid*1e3, emf)
ax[0, 0].set_xlim([0, np.max(Macid*1e3)])
ax[0, 0].set_ylim([np.min(emf)-10, np.max(emf)+10])
ax[0, 0].set_xlabel('Acid mass / g')
ax[0, 0].set_ylabel('EMF / mV')

ax[1, 0].scatter(Macid*1e3, f1)