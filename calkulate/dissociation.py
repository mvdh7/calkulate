# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019  Matthew Paul Humphreys  (GNU GPLv3)
"""Evaluate stoichiometric equilibrium constants."""
from numpy import exp, log, sqrt

def eqConstants(tempK, pSal, concTotals):
    """Assemble a dict of dissociation constants on the Free pH scale."""
    eqConstantsT = {} # for Total scale dissociation constants
    eqConstants = {} # for Free scale dissociation constants
    eqConstantsT['C1'], eqConstantsT['C2'] = ksH2CO3_T_LDK00(tempK, pSal)
    eqConstantsT['B'] = kBOH3_T_D90a(tempK, pSal)
    eqConstantsT['w'] = kH2O_T_DSC07(tempK, pSal)
    eqConstants['S'] = kHSO4_F_D90b(tempK, pSal)
    eqConstantsT['F'] = kHF_T_PF87(tempK, pSal)
    eqConstantsT['P1'], eqConstantsT['P2'], eqConstantsT['P3'] = \
        ksH3PO4_T_DSC07(tempK, pSal)
    eqConstantsT['Si'] = kSiOH4_T_M95(tempK, pSal)
    # Get pH scale conversion factor(s) and convert all to Free pH scale:
    T2F = 1/(1 + concTotals['S']/eqConstants['S'])
    eqConstants.update({k: v*T2F for k, v in eqConstantsT.items()})
    return eqConstants

def ionicStrength(pSal):
    """Estimate ionic strength from practical salinity."""
    return 19.924*pSal/(1000 - 1.005*pSal)

def ksH2CO3_T_LDK00(tempK, pSal):
    """Carbonic acid dissociation constants [LDK00]."""
    # 2 < T < 35 degC; 19 < S < 43
    pK_T_C1 = (3633.86/tempK - 61.2172 + 9.6777*log(tempK)
        - 0.011555*pSal + 0.0001152*pSal**2) # LDK00 Eq. (16)
    pK_T_C2 = (471.78/tempK + 25.929 - 3.16967*log(tempK)
        - 0.01781*pSal + 0.0001122*pSal**2) # LDK00 Eq. (17)
    return 10**-pK_T_C1, 10**-pK_T_C2

def kBOH3_T_D90a(tempK, pSal):
    """Boric acid dissociation constant [D90a]."""
    # 0 < T < 45 degC; 5 < S < 45
    lnkB = ((-8966.90 - 2890.53*sqrt(pSal) - 77.942*pSal
        + 1.728*pSal**1.5 - 0.0996*pSal**2) / tempK
        + 148.0248 + 137.1942*sqrt(pSal) + 1.62142*pSal
        - (24.4344 + 25.085*sqrt(pSal) + 0.2474*pSal) * log(tempK)
        + 0.053105*sqrt(pSal)*tempK) # D90a Eq. (23)
    return exp(lnkB)

def kH2O_T_DSC07(tempK, pSal):
    """Water dissociation constant [DSC07]."""
    lnkH2O = (148.9652 - 13847.26/tempK - 23.6521*log(tempK)
        + (118.67/tempK - 5.977 + 1.0495*log(tempK)) * sqrt(pSal)
        - 0.01615*pSal)
    return exp(lnkH2O)

def kHSO4_F_D90b(tempK, pSal):
    """Bisulfate dissociation constant [D90b]."""
    # 0 < T < 45 degC; 5 < S < 45
    I = ionicStrength(pSal)
    lnkHSO4 = (-4276.1/tempK + 141.328 - 23.093*log(tempK)
        + (-13856/tempK + 324.57 - 47.986*log(tempK)) * sqrt(I)
        + (35474/tempK - 771.54 + 114.723*log(tempK)) * I
        - (2698/tempK) * I**1.5 + (1776/tempK) * I**2
        + log(1 - 0.001005*pSal)) # D90b Eqs. (22) & (23)
    return exp(lnkHSO4)

def kHSO4_F_WM13(tempK, pSal):
    """Bisulfate dissociation constant [WM13]."""
    # WM13 equation 29, constants from Corrigendum to WM13, Table 6
    c1 = 4.24666      
    c2 = -0.152671      
    c3 = 2.67059e-2
    c4 = -4.2128e-5
    c5 = 0.2542181
    c6 = -5.09534e-3
    c7 = 7.1589e-4
    c8 = -2.91179e-3
    c9 = 2.09968e-5
    c10 = -4.03724e-5
    log10_KK = ((c1 + c2*tempK + c3*tempK*log(tempK) + c4*tempK**2)*pSal**0.5 +
        (c5 + c6*tempK + c7*tempK*log(tempK))*pSal + (c8 + c9*tempK)*pSal**1.5
        + c10*pSal**2)
    # WM13 equation 30, constants taken directly from CRP94
    a1 = 562.69486
    a2 = -102.5154
    a3 = -1.117033e-4
    a4 = 0.2477538
    a5 = -13273.75
    log10_Ko = a1 + a2*log(tempK) + a3*tempK**2 + a4*tempK + a5/tempK
    # Dissociation constant
    pkHSO4 = -(log10_KK + log10_Ko)
    return 10.0**-pkHSO4

def kHF_T_PF87(tempK, pSal):
    """Hydrogen fluoride dissociation constant [PF87]."""
    # 9 < T < 33 degC; 10 < S < 40
    lnkHF = -(-874/tempK - 0.111*sqrt(pSal) + 9.68)
    return exp(lnkHF)

def kHF_F_DR79(tempK, pSal):
    """Hydrogen fluoride dissociation constant [DR79]."""
    # 5 < T < 35 degC; 10 < S < 48
    I = ionicStrength(pSal)
    lnkF = 1590.2/tempK - 12.641 + 1.525*sqrt(I) + log(1 - 0.001005*pSal)
    return exp(lnkF)

def ksH3PO4_T_DSC07(tempK, pSal):
    """Phosphoric acid dissociation constants [DSC07]."""
    lnkP1 = (-4576.752/tempK + 115.525 - 18.453*log(tempK)
        + (-106.736/tempK + 0.69171)*sqrt(pSal)
        + (-0.65643/tempK - 0.01844)*pSal)
    lnkP2 = (-8814.715/tempK + 172.0883 - 27.927*log(tempK)
        + (-160.34/tempK + 1.3566)*sqrt(pSal)
        + (0.37335/tempK - 0.05778)*pSal)
    lnkP3 = (-3070.75/tempK - 18.141
        + (17.27039/tempK + 2.81197)*sqrt(pSal)
        + (-44.99486/tempK - 0.09984)*pSal)
    return exp(lnkP1), exp(lnkP2), exp(lnkP3)

def kSiOH4_T_M95(tempK, pSal):
    """Silicic acid dissociation constant [M95]."""
    I = ionicStrength(pSal)
    lnkSi = (-8904.2/tempK + 117.385 - 19.334*log(tempK)
        + (-458.79/tempK + 3.5913)*sqrt(I)
        + (188.74/tempK - 1.5998)*I + (-12.1652/tempK + 0.07871)*I**2
        + log(1 - 0.001005*pSal))
    return exp(lnkSi)

def kNH3_U_BJJL08(tempK, pSal):
    """Ammonium dissociation constant [BJJL08]."""
    pKNH4 = 10.0423 - 0.0315536*tempK + 0.003071*pSal # BJJL08 Eq. (3)
    return 10**-pKNH4

def k2AMP_S_BE86(tempK, pSal):
    """2-aminopyridine dissociation constant [BE86]."""
    # 5 < T < 40 degC; 30 < S < 40
    pK_S_2AMP = (2498.31/tempK - 15.3274 + 2.4050*log(tempK)
        + (0.012929 - 2.9417e-5*tempK)*pSal) # BE86 Eq. (10)
    return 10**-pK_S_2AMP
