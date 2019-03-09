# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019  Matthew Paul Humphreys  (GNU GPLv3)

from .const import psal2Cl, RMM_B, RMM_F

def XT(psal, CT=0, PT=0, SiT=0):
    return [None, CT, BT(psal), ST(psal), FT(psal), PT, SiT]

def BT(psal):
# Estimate total borate from practical salinity in mol/kg-sw [LKB10]
    return psal * 0.1336e-3 / RMM_B

def FT(psal):
# Estimate total fluoride from practical salinity [W71] in mol/kg-sw
    return psal * 6.75e-5 / (RMM_F * psal2Cl)

# Estimate total sulfate from practical salinity [????] in mol/kg-sw
def ST(psal):
    return (0.14 / 96.061) * (psal / psal2Cl)
