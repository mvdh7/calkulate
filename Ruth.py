import calkulate as calk
import os
import numpy as np

name = 'datfiles/Ruth/'
filelist = os.listdir(name)
#test = calk.gettit.VINDTA('datfiles/Ruth/Manual_Seawater.dat')

Vsamp = 100.
Cacid = 0.09
S = 33.
CT = 2000e-6
PT = 0.1e-6
SiT = 10e-6

AT   = np.full_like(filelist,np.nan, dtype='float64')
EMF0 = np.full_like(filelist,np.nan, dtype='float64')

for ix, file in enumerate(filelist):

    datfile = name + file
    print(datfile)

    AT[ix], EMF0[ix] = calk.VINDTA.MPH(datfile,Vsamp,Cacid,S,CT,PT,SiT)['x']
    
AT *= 1e6
