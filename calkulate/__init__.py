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
