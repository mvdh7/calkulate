# Calkulate: seawater total alkalinity from titration data.
# Copyright (C) 2019  Matthew Paul Humphreys  (GNU GPLv3)
"""Solve titration data for alkalinity."""
from numpy import (exp, full, isnan, log, log10, logical_and, mean, nan,
    nanmean, size, std, zeros)
from numpy import abs as np_abs
from numpy import max as np_max
from scipy.optimize import least_squares as olsq
from scipy.stats import linregress
from . import simulate
from .constants import F, R

#====== EMF to [H+] CONVERSIONS ===============================================
def emf2h(emf, emf0, tempK):
    """Convert EMF to [H+]."""
    # DAA03 Eq. (13) with typo corrected (i.e. EMF and EMF0 switched)
    return exp((emf - emf0)*F/(R*tempK))

def h2emf(h, emf0, tempK):
    """Convert [H+] to EMF."""
    return emf0 + log(h)*R*tempK/F

def f2demf0(tempK, f):
    return log(f)*R*tempK/F

#====== GRAN ESTIMATOR FUNCTIONS ==============================================
def f1(massAcid, emf, tempK, massSample):
    """Simple Gran plot estimator function, DAA03 Eq. (10)."""
    return (massSample + massAcid)*exp(emf*F/(R*tempK))

def granGuessAlk(massAcid, f1, massSample, concAcid):
    """Simple Gran plot first guess of alkalinity."""
    grad, intY, _, _, _ = linregress(massAcid, f1)
    intX = -intY/grad
    alkGuess = intX*concAcid/massSample
    return alkGuess

def granEmf0(massAcid, emf, tempK, massSample, concAcid, alk, HSO4=0, HF=0):
    """DAA03 Eq. (11)."""
    return (emf - (R*tempK/F)*log(((massAcid*concAcid - massSample*alk)
        - massSample*(HF + HSO4))/(massSample + massAcid)))

def guessGran(massAcid, emf, tempK, massSample, concAcid):
    """Simple Gran plot first guesses for AT and EMF0."""
    f1Guess = f1(massAcid, emf, tempK, massSample)
    LGuess = logical_and(
        f1Guess > 0.1*np_max(f1Guess),
        f1Guess < 0.9*np_max(f1Guess),
    )
    alkGuess = granGuessAlk(massAcid[LGuess], f1Guess[LGuess], massSample,
        concAcid)
    emf0Guess = mean(granEmf0(massAcid[LGuess], emf[LGuess], tempK[LGuess],
        massSample, concAcid, alkGuess))
    hGuess = emf2h(emf, emf0Guess, tempK)
    pHGuess = -log10(hGuess)
    return alkGuess, emf0Guess, hGuess, pHGuess

def mu(massAcid, massSample):
    """Acid dilution factor."""
    return massSample/(massSample + massAcid)

#====== LEAST-SQUARES SOLVERS =================================================
#----- Complete Calculation ---------------------------------------------------
def _lsqfun_complete(massAcid, emf, tempK, massSample, concAcid, emf0, alk, XT,
        KX):
    xmu = mu(massAcid, massSample)
    H = emf2h(emf, emf0, tempK)
    return (simulate.AT(H, xmu, XT, KX)[0] - alk*xmu +
        massAcid*concAcid/(massAcid + massSample))

def complete(massAcid, emf, tempK, massSample, concAcid, XT, KXF):
    """Solve for alkalinity using the Complete Calculation method."""
    alkGuess, emf0Guess, _, pHGuess = guessGran(massAcid, emf, tempK,
        massSample, concAcid)
    L = logical_and(pHGuess > 3, pHGuess < 4)
    KXL = {k: v[L] for k, v in KXF.items()}
    alk_emf0 = olsq(lambda AT_emf0: _lsqfun_complete(massAcid[L], emf[L],
            tempK[L], massSample, concAcid, alk_emf0[1], alk_emf0[0], XT, KXL),
        [alkGuess, emf0Guess], x_scale=[1e-6, 1], method='lm')
    return alk_emf0

def complete_emf0(massAcid, emf, tempK, emf0, massSample, concAcid, XT, KX):
    """Complete Calculation of AT with known EMF0."""
    H  = emf2h(emf, emf0, tempK)
    pH = -log10(H)
    L = logical_and(pH > 3, pH < 4)
    xmu = mu(massAcid, massSample)
    alk = (simulate.AT(H, xmu, XT, *KX)[0] + massAcid*concAcid/(massAcid +
        massSample))/xmu
    return mean(alk[L]), std(alk[L])

#----- Dickson et al. (2003) method -------------------------------------------
def _lsqfun_DAA03(massAcid, H, massSample, concAcid, f, AT, XT, KXF):
    """DAA03 Eqs. (14) and (15)."""
    Z = 1 + XT['S']/KXF['S']
    return (AT + XT['S']/(1 + KXF['S']*Z/(f*H))
        + XT['F']/(1 + KXF['F']/(f*H))
        + ((massSample + massAcid)/massSample)*f*H/Z -
            massAcid*concAcid/massSample)

def DAA03(massAcid, emf, tempK, massSample, concAcid, XT, KXF):
    """Solve for alkalinity following Dickson et al. (2003)."""
    alkGuess, emf0Guess, hGuess, pHGuess = guessGran(massAcid, emf, tempK,
        massSample, concAcid)
    L = logical_and(pHGuess > 3, pHGuess < 3.5)
    KXL = {k: v[L] for k, v in KXF.items()}
    return olsq(lambda alk_f:
        _lsqfun_DAA03(massAcid[L], hGuess[L], massSample, concAcid, alk_f[1],
            alk_f[0], XT, KXL),
        [alkGuess, 1], x_scale=[1e-3, 1], method='lm')

# To convert Dickson's "f" to an EMF0 value:
#def DAA03_emf0:
#    return EMF0g - log(alk_f['x'][1]) * R*mean(tempK)/F

#----- Dickson (1981) method --------------------------------------------------
def _lsqfun_Dickson1981(massAcid, H, massSample, concAcid, f, AT, XT, KXF):
    """Dickson (1981) Eq. (16)."""
    return (massSample*(AT - XT['C']*(1/(1 + f*H/KXF['C1'])
        + 1/(1 + f*H/KXF['C2'])) - XT['B']/(1 + f*H/KXF['B']))
        + (massSample + massAcid)*(f*H - KXF['w']/(f*H)) - massAcid*concAcid)

def Dickson1981(massAcid, EMF, tempK, massSample, concAcid, XT, KXF):
    """Solve for alkalinity and CT following Dickson (1981)."""
    alkGuess, EMF0g, hGuess, pHGuess = guessGran(massAcid, EMF, tempK,
        massSample, concAcid)
    L = pHGuess > 5
    KXL = {k: v[L] for k, v in KXF.items()}
    alk_CT_f = olsq(lambda alk_CT_f:
        _lsqfun_Dickson1981(massAcid[L], hGuess[L], massSample, concAcid,
            alk_CT_f[2], alk_CT_f[0], alk_CT_f[1], XT, KXL),
        [alkGuess, alkGuess*0.95, 1], x_scale=[1e-3, 1e-3, 1], method='lm')
    return alk_CT_f#, EMF0g - log(alk_CT_f['x'][2]) * R * mean(tempK) / F

#====== HALF-GRAN PLOT METHOD =================================================
def halfGran(massAcid, emf, tempK, massSample, concAcid, XT, KXF,
        pHRange=[3., 4.], suppressWarnings=False):
    """Solve for alkalinity using the half-Gran method (Humphreys, 2015)."""
    xmu = mu(massAcid, massSample)
    granReps = int(20)
    stepAlk = full(granReps, nan)
    stepEmf0 = full(granReps, nan)
    granHSO4 = zeros(size(emf))
    granHF = zeros(size(emf))
    granG = full((granReps, size(emf)), nan)
    granH = full((granReps, size(emf)), nan)
    granEmf0 = full((granReps, size(emf)), nan)
    granPH = full((granReps, size(emf)), nan)
    converged = False
    for i in range(granReps):
        if i == 0:
            granG[i] = f1(massAcid, emf, tempK, massSample)
            LG = granG[i] > 0.1*np_max(granG[i])
        else:
            LG = logical_and(granPH[i-1] > pHRange[0],
                             granPH[i-1] < pHRange[1])
        stepAlk[i] = granGuessAlk(massAcid[LG], granG[i, LG], massSample,
            concAcid)
        PPC = 5e-3 # permitted % change in AT
        if i > 2:
            if np_abs(stepAlk[i] - stepAlk[i-1])/stepAlk[i] < PPC/100:
                converged = True
                break
        granEmf0[i, LG] = granEmf0(massAcid[LG], emf[LG], tempK[LG],
            massSample, concAcid, stepAlk[i], granHSO4[LG], granHF[LG])
        stepEmf0[i] = nanmean(granEmf0[i])
        granH[i] = emf2h(emf, stepEmf0[i], tempK)
        granPH[i] = -log10(granH[i])
        granBicarb = xmu*XT['C']/(granH[i]/KXF['C1'] + 1)
        granHSO4 = xmu*XT['S']/(1 + KXF['S']/granH[i])
        granHF = xmu*XT['F']/(1 + KXF['F']/granH[i])
        granBorate = xmu*XT['B']/(1 + granH[i]/KXF['B'])
        granOH = KXF['w']/granH[i]
        granPP2 = (xmu*XT['P']*(1 - KXF['P1']*KXF['P2']/(granH[i]**2))
            / (1 + KXF['P1']/granH[i] + KXF['P2']*KXF['P3']/granH[i]**2
                + KXF['P1']*KXF['P2']*KXF['P3']/granH[i]**3))
        if i < granReps-1:
            granG[i+1] = (granH[i] + granHSO4 + granHF - granBicarb
                - granOH - granBorate + granPP2)*(massSample + massAcid)
    if converged:
        finalAlk = stepAlk[~isnan(stepAlk)][-1]
        finalEmf0 = stepEmf0[~isnan(stepEmf0)][-1]
    else:
        if not suppressWarnings:
            print('Calkulate: half-Gran plot iterations did not converge!')
        finalAlk = nan
        finalEmf0 = nan
    return finalAlk, finalEmf0
