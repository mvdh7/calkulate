# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew Paul Humphreys  (GNU GPLv3)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""Seawater total alkalinity from titration data."""

__all__ = [
    'calibrate',
    'concentrations',
    'constants',
    'convert',
    'datfile',
    'density',
    'dissociation',
    'io',
    'meta',
    'plot',
    'simulate',
    'solve',
]

from . import (
    calibrate,
    concentrations,
    constants,
    convert,
    datfile,
    density,
    dissociation,
    io,
    meta,
    plot,
    simulate,
    solve,
)

# Add alias to avoid breaking old code
vindta = datfile

__author__ = 'Matthew P. Humphreys and Ruth S. Matthews'
__version__ = meta.version

# Define classes
from numpy import isscalar, full_like, ndarray, size

class Potentiometric:
    """A potentiometric titration dataset."""

    def __init__(self, volAcid, emf, tempK):
        # Check inputs are correctly formatted
        assert type(volAcid) is ndarray, '`volAcid` must be a NumPy ndarray.'
        assert type(emf) is ndarray, '`emf` must be a NumPy ndarray.'
        assert size(volAcid) == size(emf), \
            '`volAcid` and `emf` must be the same size.'
        assert size(tempK) == size(volAcid) or size(tempK) == 1, \
            ('`tempK` must be either the same size as `volAcid` and `emf` or ' +
                'a single value.')
        # Ravel and store inputs
        self.volAcid = volAcid.ravel()
        self.emf = emf.ravel()
        if size(tempK) == 1:
            tempK = full_like(volAcid, tempK)
        self.tempK = tempK.ravel()

    def __repr__(self):
        fields = vars(self).keys()
        return ('<calkulate.Potentiometric with fields: ' +
            '{}, '*(len(fields) - 1) + '{}>').format(*fields)

    # def __str__(self):
    #     # Set what's returned by print(Potentiometric())
    #     return

    def apply_buretteCorrection(self, buretteCorrection):
        assert isscalar(buretteCorrection), \
            '`buretteCorrection must be a scalar.`'
        self.volAcid *= buretteCorrection
        self.buretteCorrection_applied = buretteCorrection
        # print('Burette correction has been applied to `volAcid`.')
        # print('You may now need to recalculate everything else!')

    def add(self, **kwargs):
        """Add titration variables."""
        # Define valid input keywords and test their values' formats
        scalars = ['concAcid', 'pSal', 'volSample']
        strings = ['solver']
        dicts = ['concTotals', 'eqConstants']
        valid_kwargs = [*scalars, *strings, *dicts]
        assert all([k in valid_kwargs for k in kwargs.keys()]), \
            ('Valid inputs are: ' + '`{}`, '*(len(valid_kwargs) - 1) +
            '`{}`.').format(*valid_kwargs)
        assert all([isscalar(v) for k, v in kwargs.items() if k in scalars]), \
            ('These inputs must be scalars: ' + '`{}`, '*(len(scalars) - 1) +
            '`{}`.').format(*scalars)
        assert all([isinstance(v, str) for k, v in kwargs.items()
                if k in strings]), \
            ('These inputs must be strings: ' + '`{}`, '*(len(strings) - 1) +
            '`{}`.').format(*strings)
        assert all([isinstance(v, dict) for k, v in kwargs.items()
                if k in dicts]), \
            ('These inputs must be dicts: ' + '`{}`, '*(len(dicts) - 1) +
            '`{}`.').format(*dicts)
        # Add acid concentration for a solver in mol/l
        if 'concAcid' in kwargs.keys():
            assert 'solver' in kwargs.keys(), \
                'A `solver` must be specified to accompany each `concAcid`.'
            if 'concAcid' not in vars(self):
                self.concAcid = {}
            self.concAcid[kwargs['solver']] = kwargs['concAcid']
        # Add practical salinity for the sample
        if 'pSal' in kwargs.keys():
            self.pSal = kwargs['pSal']
        # Add sample volume in ml
        if 'volSample' in kwargs.keys():
            self.volSample = kwargs['volSample']
        # Add dicts of concentrations and dissociation constants
        if 'concTotals' in kwargs.keys():
            self.concTotals = kwargs['concTotals']
        if 'eqConstants' in kwargs.keys():
            self.eqConstants = kwargs['eqConstants']

    def checkset(self, **kwargs):
        """Check if the inputs exist and if not, add them (if supplied) or
        throw an error.
        """
        for k, v in kwargs.items():
            if k == 'solver':
                assert v in solve.allSolvers.keys(), \
                    '`solver` must be in `calkulate.solve.allSolvers.keys()`.'
            elif v is not None:
                if k == 'concAcid':
                    self.add(**{k: v, 'solver': kwargs['solver']})
                else:
                    self.add(**{k: v})
                assert k in vars(self), 'A `{}` value is required.'.format(k)

    def vol2mass(self, pSal=None, volSample=None):
        """Convert sample and acid volumes to masses."""
        self.checkset(pSal=pSal, volSample=volSample)
        self.massAcid = convert.vol2massAcid(self.volAcid, self.tempK)
        self.massSample = convert.vol2massSample(self.volSample, self.tempK[0],
            self.pSal)

    def get_concTotals(self, pSal=None, totalCarbonate=None, totalPhosphate=None,
            totalSilicate=None, totalAmmonia=None, totalH2Sulfide=None,
            WhichKs=None, WhoseTB=None):
        """Calculate concentration totals from salinity."""
        gcl = {k: v for k, v in locals().items()
            if k not in ['self', 'pSal'] and v is not None}
        self.checkset(pSal=pSal)
        self.concTotals = concentrations.concTotals(self.pSal, **gcl)

    def get_eqConstants(self, pSal=None, concTotals=None, WhichKs=None,
            WhoseKSO4=None, WhoseKF=None):
        """Calculate equilibrium constants from temperature and salinity."""
        gcl = {k: v for k, v in locals().items()
            if k not in ['self', 'pSal', 'concTotals'] and v is not None}
        self.checkset(pSal=pSal, concTotals=concTotals)
        self.eqConstants = dissociation.eqConstants(self.tempK, self.pSal,
            self.concTotals, **gcl)

    def solve_alkalinity(self, concAcid=None, solver='complete', pSal=None,
            volSample=None, concTotals=None, eqConstants=None):
        """Solve for total alkalinity."""
        self.checkset(concAcid=concAcid, pSal=pSal, solver=solver,
            volSample=volSample, concTotals=concTotals, eqConstants=eqConstants)
        if 'massAcid' not in vars(self) or 'massSample' not in vars(self):
            self.vol2mass()
        if 'concTotals' not in vars(self):
            self.get_concTotals()
        if 'eqConstants' not in vars(self):
            self.get_eqConstants()
        if 'solved' not in vars(self):
            self.solved = {}
            self.alk = {}
            self.emf0 = {}
        self.solved[solver] = solve.allSolvers[solver](self.massAcid, self.emf,
            self.tempK, self.massSample, self.concAcid[solver], self.concTotals,
            self.eqConstants)
        self.alk[solver], self.emf0[solver] = self.solved[solver]['x']

    def get_alkSteps(self, solver='complete'):
        """Calculate total alkalinity separately at each titration step."""
        self.h = solve.emf2h(self.emf, self.emf0[solver], self.tempK)
        self.mu = solve.mu(self.massAcid, self.massSample)
        self.alkSim, self.alkSimComponents = simulate.alk(self.h, self.mu,
            self.concTotals, self.eqConstants)
        self.alkSteps = ((self.alkSim + self.massAcid*self.concAcid[solver]/
            (self.massAcid + self.massSample))/self.mu)
