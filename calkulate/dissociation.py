# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019  Matthew Paul Humphreys  (GNU GPLv3)

from numpy import exp, log, sqrt

def KX_F(tempK, psal, ST, FT):
    """Assemble a list of dissociation constants on the Free pH scale."""
    # Evaluate dissociation coefficients
    KC1_T, KC2_T = KC_T_LDK00(tempK, psal)
    KB_T = KB_T_D90a(tempK, psal)
    KH2O_T = KH2O_T_DSC07(tempK, psal)
    KHSO4_F = KHSO4_F_D90b(tempK, psal)
    KHF_T = KHF_T_PF87(tempK, psal)
    KP1_T, KP2_T, KP3_T = KP_T_DSC07(tempK, psal)
    KSi_T = KSi_T_M95(tempK, psal)
    # Get pH scale conversion factors and convert to Free pH scale
    T2F = 1 / (1 + ST/KHSO4_F)
    KC1_F  = KC1_T * T2F
    KC2_F = KC2_T * T2F
    KB_F = KB_T * T2F
    KH2O_F = KH2O_T * T2F
    KHF_F = KHF_T * T2F
    KP1_F = KP1_T * T2F
    KP2_F = KP2_T * T2F
    KP3_F = KP3_T * T2F
    KSi_F = KSi_T * T2F
    return [KC1_F, KC2_F, KB_F, KH2O_F, KHSO4_F, KHF_F, KP1_F, KP2_F, KP3_F,
        KSi_F]

def Istr(psal):
    """Estimate ionic strength from salinity."""
    return 19.924 * psal / (1000 - 1.005 * psal)

def KC_T_LDK00(tempK, psal):
    """Carbonic acid dissociation constants (Lueker et al., 2000)."""
    # 2 < T < 35 degC; 19 < S < 43
    pK_T_C1 = 3633.86/tempK - 61.2172 + 9.6777*log(tempK) \
        - 0.011555*psal + 0.0001152*psal**2 # LDK00 Eq. (16)
    pK_T_C2 =  471.78/tempK + 25.929 - 3.16967*log(tempK) \
        - 0.01781*psal + 0.0001122*psal**2 # LDK00 Eq. (17)
    return 10**-pK_T_C1, 10**-pK_T_C2

def KB_T_D90a(tempK, psal):
    """Boric acid dissociation constant (Dickson, 1990a)."""
    # 0 < T < 45 degC; 5 < S < 45
    ln_KB = (-8966.90 - 2890.53*sqrt(psal) - 77.942*psal \
        + 1.728*psal**1.5 - 0.0996*psal**2) / tempK \
        + 148.0248 + 137.1942*sqrt(psal) + 1.62142*psal \
        - (24.4344 + 25.085*sqrt(psal) + 0.2474*psal) * log(tempK) \
        + 0.053105*sqrt(psal)*tempK # D90a Eq. (23)
    return exp(ln_KB)

def KH2O_T_DSC07(tempK, psal):
    """Water dissociation constant (Dickson et al., 2007)."""
    ln_KH2O = 148.9652 - 13847.26/tempK - 23.6521*log(tempK) \
        + (118.67/tempK - 5.977 + 1.0495*log(tempK)) * sqrt(psal) \
        - 0.01615*psal
    return exp(ln_KH2O)

def KHSO4_F_D90b(tempK, psal):
    """Bisulfate dissociation constant (Dickson, 1990b)."""
    # 0 < T < 45 degC; 5 < S < 45
    I = Istr(psal)
    ln_KHSO4 = -4276.1/tempK + 141.328 - 23.093*log(tempK) \
        + (-13856/tempK + 324.57 - 47.986*log(tempK)) * sqrt(I) \
        + (35474/tempK - 771.54 + 114.723*log(tempK)) * I \
        - (2698/tempK) * I**1.5 + (1776/tempK) * I**2 \
        + log(1 - 0.001005*psal) # D90b Eqs. (22) & (23)
    return exp(ln_KHSO4)

def KHF_T_PF87(tempK, psal):
    """Hydrogen fluoride dissociation constant (Perez and Fraga, 1987)."""
    # 9 < T < 33 degC; 10 < S < 40
    ln_KHF = -(-874/tempK - 0.111*sqrt(psal) + 9.68)
    return exp(ln_KHF)

def KHF_F_DR79(tempK, psal):
    """Hydrogen fluoride dissociation constant (Dickson and Riley, 1979)."""
    # 5 < T < 35 degC; 10 < S < 48
    I = Istr(psal)
    ln_KF = 1590.2/tempK - 12.641 + 1.525*sqrt(I) + log(1 - 0.001005*psal)
    return exp(ln_KF)

def KP_T_DSC07(tempK, psal):
    """Phosphoric acid dissociation constants (Dickson et al., 2007)."""
    ln_KP1 = - 4576.752/tempK + 115.525 - 18.453*log(tempK) \
        + (-106.736/tempK + 0.69171) * sqrt(psal) \
        + (-0.65643/tempK - 0.01844) * psal
    ln_KP2 = - 8814.715/tempK + 172.0883 - 27.927*log(tempK) \
        + (-160.34/tempK + 1.3566) * sqrt(psal) \
        + (0.37335/tempK - 0.05778) * psal
    ln_KP3 = -3070.75/tempK - 18.141 \
        + (17.27039/tempK +  2.81197) * sqrt(psal) \
        + (-44.99486/tempK -  0.09984) * psal
    return exp(ln_KP1), exp(ln_KP2), exp(ln_KP3)

def KSi_T_M95(tempK, psal):
    """Silicic acid dissociation constant (Millero, 1995)."""
    I = Istr(psal)
    ln_KSi = -8904.2/tempK + 117.385 - 19.334*log(tempK) \
        + (-458.79/tempK + 3.5913) * sqrt(I) \
        + (188.74/tempK - 1.5998) * I + (-12.1652/tempK + 0.07871) * I**2 \
        + log(1 - 0.001005*psal)
    return exp(ln_KSi)

def KNH4_X_BJJL08(tempK, psal):
    """Ammonium dissociation constant (Bell et al., 2008)."""
    pKNH4 = 10.0423 - 0.0315536*tempK + 0.003071*psal # BJJL08 Eq. (3)
    return 10**-pKNH4

def K2AMP_S_BE86(tempK, psal):
    """2-aminopyridine dissociation constant (Bates and Erickson, 1986)."""
    # 5 < T < 40 degC; 30 < S < 40
    pK_S_2AMP = 2498.31/tempK - 15.3274 + 2.4050*log(tempK) \
        + (0.012929 - 2.9417e-5*tempK) * psal # BE86 Eq. (10)
    return 10**-pK_S_2AMP
