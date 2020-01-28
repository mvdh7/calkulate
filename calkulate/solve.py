# Calkulate: seawater total alkalinity from titration data.
# Copyright (C) 2019-2020  Matthew Paul Humphreys  (GNU GPLv3)
"""Solve titration data for alkalinity."""
from numpy import (exp, full, isnan, log, log10, logical_and, mean, nan,
    nanmean, size, zeros)
from numpy import abs as np_abs
from numpy import max as np_max
from scipy.optimize import least_squares as olsq
from scipy.stats import linregress
from . import simulate
from .constants import F, R

#====== EMF CONVERSIONS =======================================================
def emf2h(emf, emf0, tempK):
    """Convert EMF to [H+]."""
    # DAA03 Eq. (13) with typo corrected (i.e. EMF and EMF0 switched)
    return exp((emf - emf0)*F/(R*tempK))

def h2emf(h, emf0, tempK):
    """Convert [H+] to EMF."""
    return emf0 + log(h)*R*tempK/F

def f2dEmf0(tempK, f):
    return log(f)*R*tempK/F

#====== GRAN ESTIMATOR FUNCTIONS ==============================================
def f1(massAcid, emf, tempK, massSample):
    """Simple Gran plot estimator function, DAA03 Eq. (10)."""
    return (massSample + massAcid)*exp(emf*F/(R*tempK))

def granAlkGuess(massAcid, f1, massSample, concAcid):
    """Simple Gran plot first guess of alkalinity."""
    grad, intY, _, _, _ = linregress(massAcid, f1)
    intX = -intY/grad
    alkGuess = intX*concAcid/massSample
    return alkGuess

def granEmf0Guess(massAcid, emf, tempK, massSample, concAcid, alk, HSO4=0,
        HF=0):
    """DAA03 equation (11)."""
    return (emf - (R*tempK/F)*log(((massAcid*concAcid - massSample*alk)
        - massSample*(HF + HSO4))/(massSample + massAcid)))

def guessGran(massAcid, emf, tempK, massSample, concAcid):
    """Simple Gran plot first guesses for alkalinity and EMF0."""
    f1Guess = f1(massAcid, emf, tempK, massSample)
    LGuess = logical_and(
        f1Guess > 0.1*np_max(f1Guess),
        f1Guess < 0.9*np_max(f1Guess),
    )
    alkGuess = granAlkGuess(massAcid[LGuess], f1Guess[LGuess], massSample,
        concAcid)
    emf0Guess = mean(granEmf0Guess(massAcid[LGuess], emf[LGuess],
        tempK[LGuess], massSample, concAcid, alkGuess))
    hGuess = emf2h(emf, emf0Guess, tempK)
    pHGuess = -log10(hGuess)
    return alkGuess, emf0Guess, hGuess, pHGuess

def mu(massAcid, massSample):
    """Acid dilution factor."""
    return massSample/(massSample + massAcid)

#====== LEAST-SQUARES SOLVERS =================================================
def _eqConcL(concTotals, eqConstants, L):
    eqConstantsL = {k: v[L] for k, v in eqConstants.items()}
    concTotalsL = {k: v if size(v)==1 else v[L] for k, v in concTotals.items()}
    return concTotalsL, eqConstantsL

#----- Complete Calculation ---------------------------------------------------
def _lsqfunComplete(massAcid, emf, tempK, massSample, concAcid, emf0, alk,
        concTotals, eqConstants):
    xmu = mu(massAcid, massSample)
    h = emf2h(emf, emf0, tempK)
    return (simulate.alk(h, xmu, concTotals, eqConstants)[0] - alk*xmu +
        massAcid*concAcid/(massAcid + massSample))

def complete(massAcid, emf, tempK, massSample, concAcid, concTotals,
        eqConstants):
    """Solve for alkalinity and EMF0 using the complete calculation method."""
    alkGuess, emf0Guess, _, pHGuess = guessGran(massAcid, emf, tempK,
        massSample, concAcid)
    L = logical_and(pHGuess > 3, pHGuess < 4)
    concTotalsL, eqConstantsL = _eqConcL(concTotals, eqConstants, L)
    alk_emf0 = olsq(lambda alk_emf0: _lsqfunComplete(massAcid[L], emf[L],
            tempK[L], massSample, concAcid, alk_emf0[1], alk_emf0[0],
            concTotalsL, eqConstantsL),
        [alkGuess, emf0Guess], x_scale=[1e-6, 1], method='lm')
    return alk_emf0

#----- Dickson et al. (2003) method -------------------------------------------
def _lsqfun_DAA03(massAcid, H, massSample, concAcid, f, AT, concTotals,
        eqConstants):
    """DAA03 Eqs. (14) and (15)."""
    Z = 1 + concTotals['S']/eqConstants['S']
    return (AT + concTotals['S']/(1 + eqConstants['S']*Z/(f*H))
        + concTotals['F']/(1 + eqConstants['F']/(f*H))
        + ((massSample + massAcid)/massSample)*f*H/Z -
            massAcid*concAcid/massSample)

def DAA03(massAcid, emf, tempK, massSample, concAcid, concTotals, eqConstants):
    """Solve for alkalinity and f using the Dickson CRM method [DAA03]."""
    alkGuess, emf0Guess, hGuess, pHGuess = guessGran(massAcid, emf, tempK,
        massSample, concAcid)
    L = logical_and(pHGuess > 3, pHGuess < 3.5)
    concTotalsL, eqConstantsL = _eqConcL(concTotals, eqConstants, L)
    return olsq(lambda alk_f:
        _lsqfun_DAA03(massAcid[L], hGuess[L], massSample, concAcid, alk_f[1],
            alk_f[0], concTotalsL, eqConstantsL),
        [alkGuess, 1], x_scale=[1e-3, 1], method='lm')

#----- Dickson (1981) method --------------------------------------------------
def _lsqfun_Dickson1981(massAcid, h, massSample, concAcid, f, alk,
        totalCarbonate, concTotals, eqConstants):
    """D81 equation (16)."""
    return (massSample*(alk - totalCarbonate*(1/(1 + f*h/eqConstants['C1'])
        + 1/(1 + f*h/eqConstants['C2']))
        - concTotals['B']/(1 +f*h/eqConstants['B']))
        + (massSample + massAcid)*(f*h - eqConstants['w']/(f*h))
        - massAcid*concAcid)

def Dickson1981(massAcid, emf, tempK, massSample, concAcid, concTotals,
        eqConstants):
    """Solve for alkalinity, DIC and f using the closed-cell method [D81]."""
    alkGuess, EMF0g, hGuess, pHGuess = guessGran(massAcid, emf, tempK,
        massSample, concAcid)
    L = pHGuess > 5
    concTotalsL, eqConstantsL = _eqConcL(concTotals, eqConstants, L)
    alk_totalCarbonate_f = olsq(lambda alk_totalCarbonate_f:
        _lsqfun_Dickson1981(massAcid[L], hGuess[L], massSample, concAcid,
            alk_totalCarbonate_f[2], alk_totalCarbonate_f[0],
            alk_totalCarbonate_f[1], concTotalsL, eqConstantsL),
        [alkGuess, alkGuess*0.95, 1], x_scale=[1e-3, 1e-3, 1], method='lm')
    return alk_totalCarbonate_f

#====== HALF-GRAN PLOT METHOD =================================================
def halfGran(massAcid, emf, tempK, massSample, concAcid, concTotals,
        eqConstants, pHRange=[3., 4.], suppressWarnings=False):
    """Solve for alkalinity and EMF0 using the half-Gran method [H15]."""
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
        stepAlk[i] = granAlkGuess(massAcid[LG], granG[i, LG], massSample,
            concAcid)
        PPC = 5e-3 # permitted % change in AT
        if i > 2:
            if np_abs(stepAlk[i] - stepAlk[i-1])/stepAlk[i] < PPC/100:
                converged = True
                break
        granEmf0[i, LG] = granEmf0Guess(massAcid[LG], emf[LG], tempK[LG],
            massSample, concAcid, stepAlk[i], granHSO4[LG], granHF[LG])
        stepEmf0[i] = nanmean(granEmf0[i])
        granH[i] = emf2h(emf, stepEmf0[i], tempK)
        granPH[i] = -log10(granH[i])
        granBicarb = xmu*concTotals['C']/(granH[i]/eqConstants['C1'] + 1)
        granHSO4 = xmu*concTotals['S']/(1 + eqConstants['S']/granH[i])
        granHF = xmu*concTotals['F']/(1 + eqConstants['F']/granH[i])
        granBorate = xmu*concTotals['B']/(1 + granH[i]/eqConstants['B'])
        granOH = eqConstants['w']/granH[i]
        granPP2 = (xmu*concTotals['P']*(1 -
            eqConstants['P1']*eqConstants['P2']/(granH[i]**2))
            / (1 + eqConstants['P1']/granH[i] +
            eqConstants['P2']*eqConstants['P3']/granH[i]**2 +
            eqConstants['P1']*eqConstants['P2']*eqConstants['P3']/granH[i]**3))
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
    return {'x': [finalAlk, finalEmf0]}

# Dict of all solvers
allSolvers = {
    'complete': complete,
    'daa03': DAA03,
    'dickson1981': Dickson1981,
    'halfgran': halfGran,
}
