# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019  Matthew Paul Humphreys  (GNU GPLv3)
"""Evaluate stoichiometric equilibrium constants."""
from numpy import exp, log, sqrt

def KXF(tempK, pSal, XT):
    """Assemble a dict of dissociation constants on the Free pH scale."""
    KXT = {} # for Total scale dissociation constants
    KXF = {} # for Free scale dissociation constants
    # Evaluate dissociation coefficients:
    KXT['C1'], KXT['C2'] = KC_T_LDK00(tempK, pSal)
    KXT['B'] = KB_T_D90a(tempK, pSal)
    KXT['w'] = KH2O_T_DSC07(tempK, pSal)
    KXF['S'] = KHSO4_F_D90b(tempK, pSal)
    KXT['F'] = KHF_T_PF87(tempK, pSal)
    KXT['P1'], KXT['P2'], KXT['P3'] = KP_T_DSC07(tempK, pSal)
    KXT['Si'] = KSi_T_M95(tempK, pSal)
    # Get pH scale conversion factor(s) and convert all to Free pH scale:
    T2F = 1 / (1 + XT['S']/KXF['S'])
    KXF.update({k: v*T2F for k, v in KXT.items()})
    return KXF

def Istr(pSal):
    """Estimate ionic strength from salinity."""
    return 19.924*pSal/(1000 - 1.005*pSal)

def KC_T_LDK00(tempK, pSal):
    """Carbonic acid dissociation constants (Lueker et al., 2000)."""
    # 2 < T < 35 degC; 19 < S < 43
    pK_T_C1 = (3633.86/tempK - 61.2172 + 9.6777*log(tempK)
        - 0.011555*pSal + 0.0001152*pSal**2) # LDK00 Eq. (16)
    pK_T_C2 = (471.78/tempK + 25.929 - 3.16967*log(tempK)
        - 0.01781*pSal + 0.0001122*pSal**2) # LDK00 Eq. (17)
    return 10**-pK_T_C1, 10**-pK_T_C2

def KB_T_D90a(tempK, pSal):
    """Boric acid dissociation constant (Dickson, 1990a)."""
    # 0 < T < 45 degC; 5 < S < 45
    ln_KB = ((-8966.90 - 2890.53*sqrt(pSal) - 77.942*pSal
        + 1.728*pSal**1.5 - 0.0996*pSal**2) / tempK
        + 148.0248 + 137.1942*sqrt(pSal) + 1.62142*pSal
        - (24.4344 + 25.085*sqrt(pSal) + 0.2474*pSal) * log(tempK)
        + 0.053105*sqrt(pSal)*tempK) # D90a Eq. (23)
    return exp(ln_KB)

def KH2O_T_DSC07(tempK, pSal):
    """Water dissociation constant (Dickson et al., 2007)."""
    ln_KH2O = (148.9652 - 13847.26/tempK - 23.6521*log(tempK)
        + (118.67/tempK - 5.977 + 1.0495*log(tempK)) * sqrt(pSal)
        - 0.01615*pSal)
    return exp(ln_KH2O)

def KHSO4_F_D90b(tempK, pSal):
    """Bisulfate dissociation constant (Dickson, 1990b)."""
    # 0 < T < 45 degC; 5 < S < 45
    I = Istr(pSal)
    ln_KHSO4 = (-4276.1/tempK + 141.328 - 23.093*log(tempK)
        + (-13856/tempK + 324.57 - 47.986*log(tempK)) * sqrt(I)
        + (35474/tempK - 771.54 + 114.723*log(tempK)) * I
        - (2698/tempK) * I**1.5 + (1776/tempK) * I**2
        + log(1 - 0.001005*pSal)) # D90b Eqs. (22) & (23)
    return exp(ln_KHSO4)

def KHF_T_PF87(tempK, pSal):
    """Hydrogen fluoride dissociation constant (Perez and Fraga, 1987)."""
    # 9 < T < 33 degC; 10 < S < 40
    ln_KHF = -(-874/tempK - 0.111*sqrt(pSal) + 9.68)
    return exp(ln_KHF)

def KHF_F_DR79(tempK, pSal):
    """Hydrogen fluoride dissociation constant (Dickson and Riley, 1979)."""
    # 5 < T < 35 degC; 10 < S < 48
    I = Istr(pSal)
    ln_KF = 1590.2/tempK - 12.641 + 1.525*sqrt(I) + log(1 - 0.001005*pSal)
    return exp(ln_KF)

def KP_T_DSC07(tempK, pSal):
    """Phosphoric acid dissociation constants (Dickson et al., 2007)."""
    ln_KP1 = (-4576.752/tempK + 115.525 - 18.453*log(tempK)
        + (-106.736/tempK + 0.69171)*sqrt(pSal)
        + (-0.65643/tempK - 0.01844)*pSal)
    ln_KP2 = (-8814.715/tempK + 172.0883 - 27.927*log(tempK)
        + (-160.34/tempK + 1.3566)*sqrt(pSal)
        + (0.37335/tempK - 0.05778)*pSal)
    ln_KP3 = (-3070.75/tempK - 18.141
        + (17.27039/tempK + 2.81197)*sqrt(pSal)
        + (-44.99486/tempK - 0.09984)*pSal)
    return exp(ln_KP1), exp(ln_KP2), exp(ln_KP3)

def KSi_T_M95(tempK, pSal):
    """Silicic acid dissociation constant (Millero, 1995)."""
    I = Istr(pSal)
    ln_KSi = (-8904.2/tempK + 117.385 - 19.334*log(tempK)
        + (-458.79/tempK + 3.5913)*sqrt(I)
        + (188.74/tempK - 1.5998)*I + (-12.1652/tempK + 0.07871)*I**2
        + log(1 - 0.001005*pSal))
    return exp(ln_KSi)

def KNH4_X_BJJL08(tempK, pSal):
    """Ammonium dissociation constant (Bell et al., 2008)."""
    pKNH4 = 10.0423 - 0.0315536*tempK + 0.003071*pSal # BJJL08 Eq. (3)
    return 10**-pKNH4

def K2AMP_S_BE86(tempK, pSal):
    """2-aminopyridine dissociation constant (Bates and Erickson, 1986)."""
    # 5 < T < 40 degC; 30 < S < 40
    pK_S_2AMP = (2498.31/tempK - 15.3274 + 2.4050*log(tempK)
        + (0.012929 - 2.9417e-5*tempK)*pSal) # BE86 Eq. (10)
    return 10**-pK_S_2AMP
