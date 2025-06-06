# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019--2025  Matthew P. Humphreys  (GNU GPLv3)
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

from .classes import Dataset, Titration, to_Titration
from .dataset import calibrate, solve
from .meta import __author__, __version__, hello
from .read.metadata import (
    read_clipboard,
    read_csv,
    read_dbs,
    read_excel,
    read_fwf,
    read_table,
)
from .read.titrations import read_dat, write_dat


# For backwards-compatibility
say_hello = hello
__all__ = [
    "__author__",
    "__version__",
    "hello",
    "say_hello",
    "Dataset",
    "calibrate",
    "solve",
    "read_dat",
    "read_clipboard",
    "read_csv",
    "read_dbs",
    "read_excel",
    "read_fwf",
    "read_table",
    "write_dat",
    "Titration",
    "to_Titration",
]
