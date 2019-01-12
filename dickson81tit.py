import numpy as np
import pandas as pd

tdata = pd.read_csv('dickson81tit.csv')

mass_sample = np.float_(200e-3) # kg
conc_acid = np.float_(0.3) # mol/kg-sol'n

TA = np.float_(0.00245) # Alkalinity / mol/kg-sw
TB = np.float_(0.00042) # Borate     / mol/kg-sw
TC = np.float_(0.00220) # Carbon     / mol/kg-sw
TS = np.float_(0.02824) # Sulfate    / mol/kg-sw
TF = np.float_(0.00007) # Fluoride   / mol/kg-sw
TP = np.float_(0      ) # Phosphate  / mol/kg-sw

# Constants on Free H scale
Kw    = np.float_(4.32e-14)
K1    = np.float_(1.00e-06)
K1K2  = np.float_(8.20e-16)
K2    = K1K2 / K1
KB    = np.float_(1.78e-09)
bHSO4 = np.float_(1.23e+01)
bHF   = np.float_(4.08e+02)

# Physical properties
sal = np.float_(35)
tmp = np.float_(25)
tmk = tmp + np.float_(273.15)
