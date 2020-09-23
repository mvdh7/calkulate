# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew P. Humphreys  (GNU GPLv3)
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
"""Determine (sea)water total alkalinity from titration data."""

from . import (
    constants,
    convert,
    datasets,
    density,
    options,
    simulate,
    solvers,
    titrations,
)
from .titrations import Titration, read_dat, to_dat, kwargs_TiTouch
from .datasets import (
    Dataset,
    read_csv,
    read_dbs,
    read_excel,
    get_titrations,
    get_analyte_temperature,
    get_analyte_mass,
    get_analyte_totals,
    get_titration_totals,
    get_totals,
    get_k_constants,
    set_batch_mean_molinity,
    prepare,
    solve,
    solve_all,
    calibrate,
    calibrate_all,
    calkulate,
)

calibrate_and_solve = calkulate


# Package metadata
_authorlist = ["Humphreys, Matthew P.", "Matthews, Ruth S."]
__author__ = " and ".join(_authorlist)
__version__ = "3.0.1"


def say_hello():
    """Report the version number."""
    print(r"""
   .--.     . .         .      .      
  :         | |         |     _|_     
  |    .-.  | |.-. .  . | .-.  |  .-. 
  :   (   ) | |-.' |  | |(   ) | (.-' 
   `--'`-'`-`-'  `-`--`-`-`-'`-`-'`--'

       doi:10.5281/zenodo.2634304
       
             Version {}
""".format(__version__))
