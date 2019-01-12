import calkulate as calk
import numpy as np
import os
import regex as re

datfiles = os.listdir('datfiles')

reCRM = re.compile('0-0  0  \(0\)CRM-(\d{3})-\d{4}-\d\.dat')

DicksonCRM = {144: np.array([2031.53 ,
                             2238.60 , 
                               33.571, 
                               0.31  , 
                               2.5   , 
                               0.00  ,  
                               0.85  ]),
              151: np.array([2033.83 ,
                             2225.56 ,
                               33.345,
                                0.56 ,
                                3.5  ,
                                0.00 ,
                                1.53 ]),
              160: np.array([2030.39 ,
                             2212.44 ,
                               33.414,
                                0.57 ,
                                8.2  ,
                                0.01 ,
                                1.50 ]),
              161: np.array([2037.38 ,
                             2207.33 ,
                               33.356,
                                0.59 ,
                                4.3  ,
                                0.00 ,
                                3.30 ])}

CRMbatch = np.full(len(datfiles),np.nan)
S        = np.full(len(datfiles),np.nan)
CT       = np.full(len(datfiles),np.nan)
PT       = np.full(len(datfiles),np.nan)
SiT      = np.full(len(datfiles),np.nan)
AT_cert  = np.full(len(datfiles),np.nan)

for i,datfile in enumerate(datfiles):
    
    CRMbatch[i] = int(reCRM.match(datfile).group(1))
    CT[i],AT_cert[i],S[i],PT[i],SiT[i],_,_ = DicksonCRM[CRMbatch[i]]

CT      *= 1e-6
AT_cert *= 1e-6
PT      *= 1e-6
SiT     *= 1e-6

Vsamp = 100 # ml

#%% Calibrate
Cacids = np.full(len(datfiles),np.nan)
EMF0s  = np.full(len(datfiles),np.nan)

for i,datfile in enumerate(datfiles):
    
    Cacids[i],_,EMF0s[i] = calk.VINDTA.MPH_CRM('datfiles/' + datfile,
        Vsamp,AT_cert[i],S[i],CT[i],PT[i],SiT[i])
   
# Get best calibrated values
Cacid = np.median(Cacids)
EMF0  = np.median(EMF0s )
    
# Solve for AT
AT_varE0 = np.full(len(datfiles),np.nan)
E0_varE0 = np.full(len(datfiles),np.nan)

for i,datfile in enumerate(datfiles):
    
    AT_varE0[i],E0_varE0[i] = calk.VINDTA.MPH('datfiles/' + datfile,
        Vsamp,Cacid,S[i],CT[i],PT[i],SiT[i])['x']

AT_varE0_diff = AT_varE0 - AT_cert

#%% Test fixed E0

AT_fixE0     = np.full(len(datfiles),np.nan)
AT_fixE0_std = np.full(len(datfiles),np.nan)

for i,datfile in enumerate(datfiles):

    Macid,EMF,Tk,Msamp,XT,KX = calk.VINDTA.prep('datfiles/' + datfile,
                                                Vsamp,S[i],CT[i],PT[i],SiT[i])

    AT_fixE0[i],AT_fixE0_std[i] = calk.solve.MPH_E0(Macid,EMF,Tk,EMF0,
                                                    Msamp,Cacid,XT,KX)
    
AT_fixE0_diff = AT_fixE0 - AT_cert

#%% Show results
print('Variable E0: std = %5.2f' % np.std(AT_varE0_diff*1e6))
print('Fixed    E0: std = %5.2f' % np.std(AT_fixE0_diff*1e6))
