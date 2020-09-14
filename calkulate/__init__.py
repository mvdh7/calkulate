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

from . import (
    constants,
    convert,
    datasets,
    density,
    options,
    simulate,
    solve,
    titrations,
)
from .titrations import Titration, read_dat, to_dat
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
    prepare,
)


# Package metadata
__author__ = "Matthew P. Humphreys"
__version__ = "3.0.0-beta.3"


def say_hello():
    """Report the version number."""
    print("This is Calkulate v{}.".format(__version__))
