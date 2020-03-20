# Calkulate: seawater total alkalinity from titration data.
# Copyright (C) 2019-2020  Matthew Paul Humphreys  (GNU GPLv3)
"""Function wrappers for working with titration data in text files.

# Common inputs:

  * `datFile` = .dat file name (and path).
  * `volSample` = sample volume, ml.
  * `concAcid` = acid molality, mol/kg.
  * `pSal` = practical salinity, dimensionless.
  * `alkCert` = certified total alkalinity, mol/kg-sw.
  * `totalCarbonate` = dissolved inorganic carbon, mol/kg-sw.
  * `totalPhosphate` = total phosphate, mol/kg-sw.
  * `totalSilicate` = total silicate, mol/kg-sw.
  * `tempKForce` (optional) = titration temperature, K.
  * `totalAmmonia` (optional, default 0) = total ammonia, mol/kg-sw.
  * `totalH2Sulfide` (optional, default 0) = total hydrogen sulfide, mol/kg-sw.

The final four inputs (all optional) determine how the equilibrium constants
and concentrations estimated from salinity are generated. These are calculated
by PyCO2SYS (see https://github.com/mvdh7/PyCO2SYS). These first two are
identical to PyCO2SYS inputs:

  * `WhichKs` (optional, default 10) = `K1K2CONSTANTS` in PyCO2SYS.
  * `WhoseKF` (optional, default 1) = `KFCONSTANT` in PyCO2SYS.

The other two split up PyCO2SYS's `KSO4CONSTANTS` input, which controls both
the bisulfate dissociation constant and the total borate:chlorinity ratio:

  * `WhoseKSO4` (optional, default 1) = bisulfate dissociation constant.
  * `WhoseTB` (optional, default 2) = total borate:chlorinity ratio.

The equivalent options are:

|   PyCO2SYS    | Calkulate | Calkulate |
| KSO4CONSTANTS | WhoseKSO4 |  WhoseTB  |
|:-------------:|:---------:|:---------:|
|      1        |     1     |     1     |
|      2        |     2     |     1     |
|      3        |     1     |     2     |
|      4        |     2     |     2     |
"""
from . import calibrate, concentrations, density, dissociation, io, solve

def prep(datFile, volSample, pSal, totalCarbonate, totalPhosphate,
        totalSilicate, buretteCorrection=1, tempKForce=None, WhichKs=10,
        WhoseKSO4=1, WhoseKF=1, WhoseTB=2, totalAmmonia=0, totalH2Sulfide=0):
    """Import .dat file and prepare data for analysis."""
    volAcid, emf, tempK = io.datfile(datFile)
    if tempKForce is not None:
        tempK[:] = tempKForce
    massSample = volSample*density.sw(tempK[0], pSal)*1e-3
    massAcid = buretteCorrection*volAcid*density.acid(tempK)*1e-3
    concTotals = concentrations.concTotals(pSal, totalCarbonate=totalCarbonate,
        totalPhosphate=totalPhosphate, totalSilicate=totalSilicate,
        WhichKs=WhichKs, WhoseTB=WhoseTB, totalAmmonia=totalAmmonia,
        totalH2Sulfide=totalH2Sulfide)
    eqConstants = dissociation.eqConstants(tempK, pSal, concTotals,
        WhichKs=WhichKs, WhoseKSO4=WhoseKSO4, WhoseKF=WhoseKF)
    return massAcid, emf, tempK, massSample, concTotals, eqConstants

def alk(datFile, volSample, concAcid, pSal, totalCarbonate, totalPhosphate,
        totalSilicate, solver='complete', buretteCorrection=1,
        tempKForce=None, WhichKs=10, WhoseKSO4=1, WhoseKF=1, WhoseTB=2,
        totalAmmonia=0, totalH2Sulfide=0, **kwargs):
    """Solve for alkalinity from a titration .dat file."""
    massAcid, emf, tempK, massSample, concTotals, eqConstants = prep(datFile,
        volSample, pSal, totalCarbonate, totalPhosphate, totalSilicate,
        buretteCorrection=buretteCorrection, tempKForce=tempKForce,
        WhichKs=WhichKs, WhoseKSO4=WhoseKSO4, WhoseKF=WhoseKF, WhoseTB=WhoseTB,
        totalAmmonia=totalAmmonia, totalH2Sulfide=totalH2Sulfide)
    if solver.lower() in solve.allSolvers.keys():
        solveFunc = solve.allSolvers[solver.lower()]
        alkOptResult = solveFunc(massAcid, emf, tempK, massSample, concAcid,
            concTotals, eqConstants, **kwargs)
    else:
        print('calkulate.datfile.alk: solver not recognised.')
        print('Options (case-insensitive):' +
            (len(solve.allSolvers.keys())*' \'{}\'').format(
                *solve.allSolvers.keys()))
        alkOptResult = {'x': [None,]}
    return alkOptResult

def concAcid(datFile, volSample, alkCert, pSal, totalCarbonate,
        totalPhosphate, totalSilicate, solver='complete', buretteCorrection=1,
        tempKForce=None, WhichKs=10, WhoseKSO4=1, WhoseKF=1, WhoseTB=2,
        totalAmmonia=0, totalH2Sulfide=0, **kwargs):
    """Solve for acid concentration from a titration .dat file."""
    massAcid, emf, tempK, massSample, concTotals, eqConstants = prep(datFile,
        volSample, pSal, totalCarbonate, totalPhosphate, totalSilicate,
        buretteCorrection=buretteCorrection, tempKForce=tempKForce,
        WhichKs=WhichKs, WhoseKSO4=WhoseKSO4, WhoseKF=WhoseKF, WhoseTB=WhoseTB,
        totalAmmonia=totalAmmonia, totalH2Sulfide=totalH2Sulfide)
    concAcidOptResult = calibrate.concAcid(massAcid, emf, tempK, massSample,
        alkCert, concTotals, eqConstants, solver=solver, **kwargs)
    return concAcidOptResult
