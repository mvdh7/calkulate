from .const import S_Cl, RMM_B, RMM_F

def XT(S, CT=0, PT=0, SiT=0):
    return [None, CT, BT(S), ST(S), FT(S), PT, SiT]

def BT(S):
# Estimate total boron from practical salinity in mol/kg-sw [LKB10]
    return S * 0.1336e-3 / RMM_B

def FT(S):
# Estimate total fluoride from practical salinity [W71] in mol/kg-sw
    return S * 6.75e-5 / (RMM_F * S_Cl)

# Estimate total sulfate from practical salinity [????] in mol/kg-sw
def ST(S):
    return (0.14 / 96.061) * (S / S_Cl)
