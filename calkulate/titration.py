# Calkulate: seawater total alkalinity from titration data.
# Copyright (C) 2019-2020  Matthew Paul Humphreys  (GNU GPLv3)
"""Classes for object-oriented style calculations."""
from copy import deepcopy
from numpy import isscalar, full_like, mean, nan, ndarray, size, sqrt
from numpy import all as np_all
from . import calibrate, concentrations, convert, dissociation, simulate, solve

class Potentiometric:
    """A potentiometric titration dataset."""

    def __init__(self, volAcid, emf, tempK, pSal, volSample, totalCarbonate=0,
            totalPhosphate=0, totalSilicate=0, totalAmmonia=0, totalH2Sulfide=0,
            WhichKs=10, WhoseKSO4=1, WhoseKF=1, WhoseTB=2, warnings=True,
            checkInputs=True, buretteCorrection=1):
        # Check inputs are correctly formatted, if not requested otherwise
        if checkInputs:
            assert type(volAcid) is ndarray, '`volAcid` must be an ndarray.'
            assert type(emf) is ndarray, '`emf` must be an ndarray.'
            assert size(volAcid) == size(emf), \
                '`volAcid` and `emf` must be the same size.'
            assert ((type(tempK) is ndarray and size(tempK) == size(volAcid))
                or isscalar(tempK)), \
                ('`tempK` must be either an ndarray the same size as ' +
                 '`volAcid` and `emf` or a scalar.')
            assert np_all(tempK >= 0), 'All `tempK` values must be positive.'
            assert isscalar(pSal) and pSal >= 0, \
                '`pSal` must be a scalar with a positive value.'
            assert isscalar(volSample) and volSample > 0, \
                '`volSample` must be a scalar with a positive value.'
            assert ((type(totalCarbonate) is ndarray and
                    size(totalCarbonate) == size(volAcid)) or
                    isscalar(totalCarbonate)), \
                ('`totalCarbonate` must be either an ndarray the same size ' +
                 'as `volAcid` and `emf` or a scalar.')
            assert np_all(totalCarbonate >= 0), \
                'All `totalCarbonate` values must be positive.'
            assert isscalar(totalPhosphate) and totalPhosphate >= 0, \
                '`totalPhosphate` must be a scalar with a positive value.'
            assert isscalar(totalSilicate) and totalSilicate >= 0, \
                '`totalSilicate` must be a scalar with a positive value.'
            assert isscalar(totalAmmonia) and totalAmmonia >= 0, \
                '`totalAmmonia` must be a scalar with a positive value.'
            assert isscalar(totalH2Sulfide) and totalH2Sulfide >= 0, \
                '`totalH2Sulfide` must be a scalar with a positive value.'
            assert isinstance(WhichKs, int) and 1 <= WhichKs <= 15, \
                '`WhichKs` must be a scalar integer from 1 to 15.'
            assert isinstance(WhoseKSO4, int) and 1 <= WhoseKSO4 <= 2, \
                '`WhoseKSO4` must be a scalar integer from 1 to 2.'
            assert isinstance(WhoseKF, int) and 1 <= WhoseKF <= 2, \
                '`WhoseKF` must be a scalar integer from 1 to 2.'
            assert isinstance(WhoseTB, int) and 1 <= WhoseTB <= 2, \
                '`WhoseTB` must be a scalar integer from 1 to 2.'
            assert isinstance(warnings, bool), \
                '`warnings` must be `True` or `False`.'
        # Format and store inputs
        self.volAcid = volAcid.ravel()
        self.emf = emf.ravel()
        if isscalar(tempK):
            tempK = full_like(volAcid, tempK)
        self.tempK = tempK.ravel()
        self.pSal = pSal
        self.volSample = volSample
        self.concTotals = concentrations.concTotals(pSal,
            totalCarbonate=totalCarbonate, totalPhosphate=totalPhosphate,
            totalSilicate=totalSilicate, totalAmmonia=totalAmmonia,
            totalH2Sulfide=totalH2Sulfide, WhichKs=WhichKs, WhoseTB=WhoseTB)
        self.eqConstants = dissociation.eqConstants(self.tempK, pSal,
            self.concTotals, WhichKs=WhichKs, WhoseKSO4=WhoseKSO4,
            WhoseKF=WhoseKF)
        self.WhichKs = WhichKs
        self.WhoseKSO4 = WhoseKSO4
        self.WhoseKF = WhoseKF
        self.WhoseTB = WhoseTB
        self.__warnings__ = warnings
        # Do calculations
        self.apply_buretteCorrection(buretteCorrection, warning=False)

    def __repr__(self):
        fields = vars(self).keys()
        return ('<calkulate.Potentiometric with fields: ' +
            '{}, '*(len(fields) - 1) + '{}>').format(*fields)

    def apply_buretteCorrection(self, buretteCorrection, warning=True):
        """Apply or update the burette correction factor for `volAcid`."""
        assert isscalar(buretteCorrection) and buretteCorrection > 0, \
            '`buretteCorrection` must be a positive scalar.'
        if 'buretteCorrection_applied' not in vars(self).keys():
            self.volAcid_raw = self.volAcid.copy()
        self.volAcid = self.volAcid_raw*buretteCorrection
        self.buretteCorrection_applied = buretteCorrection
        # Also apply to acid mass
        self.vol2mass()
        if self.__warnings__ and warning:
            print('Burette correction has been applied to `volAcid` and ' +
                  '`massAcid`, but not to any subsequent calculations!')

    def vol2mass(self):
        """Convert sample and acid volumes to masses."""
        self.massAcid = convert.vol2massAcid(self.volAcid, self.tempK)
        self.massSample = convert.vol2massSample(self.volSample, self.tempK[0],
            self.pSal)

    def solve_alk(self, concAcid, solver='complete', checkInputs=True):
        """Solve for total alkalinity."""
        if checkInputs:
            assert isscalar(concAcid) and concAcid > 0, \
                '`concAcid` must be a scalar with a positive value.'
            assert solver in solve.allSolvers.keys(), \
                '`solver` must be in `calkulate.solve.allSolvers.keys()`.'
        if 'alkSolved' not in vars(self):
            self.alkSolved = {}
            self.concAcid = {}
            self.alk = {}
            self.emf0 = {}
            self.solvedWith = {}
        self.concAcid[solver] = concAcid
        self.alkSolved[solver] = solve.allSolvers[solver](self.massAcid,
            self.emf, self.tempK, self.massSample, self.concAcid[solver],
            self.concTotals, self.eqConstants)
        self.alk[solver] = self.alkSolved[solver]['x'][0]
        if solver in ['complete', 'halfgran']:
            self.emf0[solver] = self.alkSolved[solver]['x'][1]
        else:
            self.emf0[solver] = nan
        self.rms = {}
        if solver in ['complete']:
            self.rms[solver] = sqrt(mean(self.alkSolved[solver]['fun']**2))
        self.solvedWith[solver] = self.alkSolved[solver]['L']

    def calibrate_concAcid(self, alkCert, solver='complete', checkInputs=True,
            **kwargs):
        """Solve for acid concentration."""
        if checkInputs:
            assert isscalar(alkCert) and alkCert > 0, \
                '`alkCert` must be a scalar with a positive value.'
            assert solver in solve.allSolvers.keys(), \
                '`solver` must be in `calkulate.solve.allSolvers.keys()`.'
        self.alkCert = alkCert
        if 'concAcid_calibrated' not in vars(self):
            self.concAcid_calibrated = {}
            self.concAcidCert = {}
        self.concAcid_calibrated[solver] = calibrate.concAcid(self.massAcid,
            self.emf, self.tempK, self.massSample, self.alkCert,
            self.concTotals, self.eqConstants, solver=solver, **kwargs)
        self.concAcidCert[solver] = self.concAcid_calibrated[solver]['x'][0]

    def get_alkSteps(self, solver='complete', checkInputs=True):
        """Calculate total alkalinity separately at each titration step."""
        if checkInputs:
            assert solver in solve.allSolvers.keys(), \
                '`solver` must be in `calkulate.solve.allSolvers.keys()`.'
        assert 'concAcid' in vars(self) and solver in self.concAcid.keys(), \
            ('You must first solve for alkalinity ' +
             'using `solve_alk()` with the same `solver`.')
        self.h = convert.emf2h(self.emf, self.emf0[solver], self.tempK)
        self.mu = solve.mu(self.massAcid, self.massSample)
        self.alkSim, self.alkComponents = simulate.alk(self.h, self.mu,
            self.concTotals, self.eqConstants)
        self.alkSteps = ((self.alkSim + self.massAcid*self.concAcid[solver]/
            (self.massAcid + self.massSample))/self.mu)

    def copy(self):
        """Return an independent copy of this titration object."""
        return deepcopy(self)
