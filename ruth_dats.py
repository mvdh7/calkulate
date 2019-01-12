import calkulate as calk
import numpy as np

Cacid = calk.VINDTA.MPH_CRM('Ruth Pre-Locate/8888_133_0709_0_0_2.dat',
                           100.,2224.37e-6,33.341,2021.12e-6,0.40e-6,2.9e-6)

crmt = calk.VINDTA.MPH('Ruth Pre-Locate/8888_133_0709_0_0_2.dat',
                       100.,Cacid[0],33.341,2021.12e-6,0.40e-6,2.9e-6)

test = calk.VINDTA.MPH('Ruth Pre-Locate/1_20190108_1_0_0_1.dat',
                       100.,Cacid[0],0.,0.,0.,0.)

ta = test['x']

Vacid, EMF, Tk = calk.gettit.VINDTA('Ruth Pre-Locate/1_20190108_1_0_0_1.dat')

EMF0 = crmt['x'][1]

H = calk.solve.EMF2H(EMF,EMF0,Tk)

pH = -np.log10(H)

from matplotlib import pyplot as plt

plt.scatter(Vacid,pH)



