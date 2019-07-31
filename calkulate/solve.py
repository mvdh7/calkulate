# Calkulate: seawater total alkalinity from titration data.
# Copyright (C) 2019  Matthew Paul Humphreys  (GNU GPLv3)
"""Solve titration datasets for alkalinity."""
from numpy import (exp, full, isnan, log, log10, logical_and, mean, nan,
    nanmean, size, std, zeros)
from numpy import abs as np_abs
from numpy import max as np_max
from scipy.optimize import least_squares as olsq
from scipy.stats import linregress
from . import simulate
from .constants import F, R

#==============================================================================
#====== EMF to [H+] CONVERSIONS ===============================================

def emf2h(emf, emf0, tempK):
    """Convert EMF to [H+]."""
    # DAA03 Eq. (13) with typo corrected (i.e. EMF and EMF0 switched)
    return exp((emf - emf0) * F/(R*tempK))

def h2emf(h, emf0, tempK):
    """Convert [H+] to EMF."""
    return emf0 + log(h) * R*tempK/F

def f2demf0(tempK, f):
    return log(f) * R*tempK/F

#==============================================================================
#====== GRAN ESTIMATOR FUNCTIONS ==============================================

def f1(Macid, emf, tempk, Msamp):
    return (Msamp + Macid) * exp(emf * F / (R * tempk)) # DAA03 Eq. (10)

def Gran_guess_AT(Macid, f1, Msamp, Cacid):
    grad, inty, _, _, _ = linregress(Macid, f1)
    intx = -inty / grad
    AT0 = intx * Cacid / Msamp
    return AT0

def Gran_EMF0(Macid, EMF, tempK, Msamp, Cacid, AT, HSO4=0, HF=0):
    return EMF - (R * tempK / F) * log(((Macid * Cacid - Msamp * AT) \
        - Msamp * (HF + HSO4)) / (Msamp + Macid)) # DAA03 Eq. (11)

def guessGran(Macid, EMF, tempK, Msamp, Cacid):
    """Get first guesses for AT and EMF0."""
    f1g = f1(Macid, EMF, tempK, Msamp)
    Lg = logical_and(f1g > 0.1 * np_max(f1g), f1g < 0.9 * np_max(f1g))
    ATg = Gran_guess_AT(Macid[Lg], f1g[Lg], Msamp, Cacid)
    EMF0g = mean(Gran_EMF0(Macid[Lg], EMF[Lg], tempK[Lg], Msamp, Cacid, ATg))
    Hg = emf2h(EMF, EMF0g, tempK)
    pHg = -log10(Hg)
    return ATg, EMF0g, Hg, pHg

#==============================================================================
#====== LEAST-SQUARES SOLVERS =================================================

#----- Complete Calculation ---------------------------------------------------

def _lsqfun_complete(Macid, EMF, tempK, Msamp, Cacid, EMF0, AT, XT, KX):
    mu = Msamp / (Msamp + Macid)
    H = emf2h(EMF, EMF0, tempK)
    return simulate.AT(H, mu, XT, KX)[0] - AT*mu + Macid*Cacid/(Macid + Msamp)

def complete(Macid, emf, tempK, Msamp, Cacid, XT, KXF):
    """Solve for AT using the Complete Calculation method."""
    ATg, emf0g, _, pHg = guessGran(Macid, emf, tempK, Msamp, Cacid)
    L = logical_and(pHg > 3, pHg < 4)
    KXL = {k: v[L] for k, v in KXF.items()}
    AT_emf0 = olsq(lambda AT_emf0: _lsqfun_complete(Macid[L], emf[L], tempK[L],
            Msamp, Cacid, AT_emf0[1], AT_emf0[0], XT, KXL),
        [ATg, emf0g], x_scale=[1e-6, 1], method='lm')
    return AT_emf0

def complete_emf0(Macid, emf, tempK, emf0, Msamp, Cacid, XT, KX):
    """Complete Calculation of AT with known EMF0."""
    H  = emf2h(emf, emf0, tempK)
    pH = -log10(H)
    L = logical_and(pH > 3, pH < 4)
    mu = Msamp / (Msamp + Macid)
    AT = (simulate.AT(H, mu, XT, *KX)[0] + Macid*Cacid/(Macid + Msamp)) / mu
    return mean(AT[L]), std(AT[L])


#----- Dickson et al. (2003) method -------------------------------------------

def _lsqfun_DAA03(Macid, H, Msamp, Cacid, f, AT, XT, KXF):
    Z = 1 + XT['S']/KXF['S'] # DAA03 Eq. (15)
    return AT + XT['S'] / (1 + KXF['S']*Z/(f*H)) \
        + XT['F']/(1 + KXF['F']/(f*H)) \
        + ((Msamp + Macid)/Msamp) * f*H/Z - Macid*Cacid/Msamp # DAA03 Eq. (14)

def DAA03(Macid, emf, tempK, Msamp, Cacid, XT, KXF):
    """Solve for AT following Dickson et al. (2003)."""    
    ATg, emf0g, Hg, pHg = guessGran(Macid, emf, tempK, Msamp, Cacid)
    L = logical_and(pHg > 3, pHg < 3.5)
    KXL = {k: v[L] for k, v in KXF.items()}
    return olsq(lambda AT_f: \
        _lsqfun_DAA03(Macid[L], Hg[L], Msamp, Cacid, AT_f[1], AT_f[0],
            XT, KXL),
        [ATg, 1], x_scale=[1e-3, 1], method='lm')

# To convert Dickson's "f" to an EMF0 value:
#def DAA03_emf0:
#    return EMF0g - log(AT_f['x'][1]) * R*mean(tempK)/F

#----- Dickson (1981) method --------------------------------------------------

def _lsqfun_Dickson1981(Macid, H, Msamp, Cacid, f, AT, XT, KXF):
    # Dickson (1981) Eq. (16):
    return Msamp * (AT - XT['C'] * (1/(1 + f*H/KXF['C1']) \
        + 1/(1 + f*H/KXF['C2'])) - XT['B']/(1 + f*H/KXF['B'])) \
        + (Msamp + Macid)*(f*H - KXF['w']/(f*H)) - Macid*Cacid

def Dickson1981(Macid, EMF, tempK, Msamp, Cacid, XT, KXF):
    """Solve for AT and CT following Dickson (1981)."""
    ATg, EMF0g, Hg, pHg = guessGran(Macid, EMF, tempK, Msamp, Cacid)
    L = pHg > 5
    KXL = {k: v[L] for k, v in KXF.items()}
    AT_CT_f = olsq(lambda AT_CT_f: \
        _lsqfun_Dickson1981(Macid[L], Hg[L], Msamp, Cacid, AT_CT_f[2],
            AT_CT_f[0], AT_CT_f[1], XT, KXL),
        [ATg,ATg*0.95, 1], x_scale=[1e-3, 1e-3, 1], method='lm')
    return AT_CT_f#, EMF0g - log(AT_CT_f['x'][2]) * R * mean(tempK) / F


#==============================================================================
#====== HALF-GRAN PLOT METHOD =================================================

def halfGran(Macid, EMF, tempK, Msamp, Cacid, XT, KXF,
        pHrange=[3., 4.], suppress_warnings=False):
    """Solve for AT using the half-Gran method (Humphreys, 2015)."""
    mu = Msamp / (Msamp + Macid)
    Gran_reps = int(20)
    step_AT = full(Gran_reps, nan)
    step_E0 = full(Gran_reps, nan)
    Gran_HSO4 = zeros(size(EMF))
    Gran_HF = zeros(size(EMF))
    Gran_G = full((Gran_reps, size(EMF)), nan)
    Gran_H = full((Gran_reps, size(EMF)), nan)
    Gran_E0 = full((Gran_reps, size(EMF)), nan)
    Gran_pH = full((Gran_reps, size(EMF)), nan)
    converged = False
    for i in range(Gran_reps):
        if i == 0:
            Gran_G[i] = f1(Macid, EMF, tempK, Msamp)
            LG = Gran_G[i] > 0.1 * np_max(Gran_G[i])
        else:
            LG = logical_and(Gran_pH[i - 1] > pHrange[0],
                             Gran_pH[i - 1] < pHrange[1])
        step_AT[i] = Gran_guess_AT(Macid[LG], Gran_G[i, LG], Msamp, Cacid)
        PPC = 5e-3 # permitted % change in AT
        if i > 2:
            if np_abs(step_AT[i] - step_AT[i - 1]) / step_AT[i] < PPC / 100:
                converged = True
                break
        Gran_E0[i, LG] = Gran_EMF0(Macid[LG], EMF[LG], tempK[LG],
            Msamp, Cacid, step_AT[i], Gran_HSO4[LG], Gran_HF[LG])
        step_E0[i] = nanmean(Gran_E0[i])
        Gran_H[i] = emf2h(EMF, step_E0[i], tempK)
        Gran_pH[i] = -log10(Gran_H[i])
        Gran_bicarb = mu * XT['C'] / (Gran_H[i]/KXF['C1'] + 1)
        Gran_HSO4 = mu * XT['S'] / (1 + KXF['S']/Gran_H[i])
        Gran_HF = mu * XT['F'] / (1 + KXF['F']/Gran_H[i])
        Gran_borate = mu * XT['B'] / (1 + Gran_H[i]/KXF['B'])
        Gran_OH = KXF['w'] / Gran_H[i]
        Gran_P_P2 = mu * XT['P'] * (1 - KXF['P1']*KXF['P2']/(Gran_H[i]**2)) \
            / (1 + KXF['P1']/Gran_H[i] + KXF['P2']*KXF['P3']/Gran_H[i]**2 \
                + KXF['P1']*KXF['P2']*KXF['P3']/Gran_H[i]**3)
        if i < Gran_reps - 1:
            Gran_G[i+1] = (Gran_H[i] + Gran_HSO4 + Gran_HF - Gran_bicarb \
                - Gran_OH - Gran_borate + Gran_P_P2) * (Msamp + Macid)
    if converged:
        final_AT = step_AT[~isnan(step_AT)][-1]
        final_E0 = step_E0[~isnan(step_E0)][-1]
    else:
        if not suppress_warnings:
            print('Calkulate: half-Gran plot iterations did not converge!')
        final_AT = nan
        final_E0 = nan
    return final_AT, final_E0
