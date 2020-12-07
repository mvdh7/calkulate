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

# Package metadata
_authorlist = ["Humphreys, Matthew P.", "Matthews, Ruth S."]
__author__ = " and ".join(_authorlist)
__version__ = "3.2.0"


def hello():
    """Report the version number."""
    print(
        r"""
   .--.     . .         .      .      
  :         | |         |     _|_     
  |    .-.  | |.-. .  . | .-.  |  .-. 
  :   (   ) | |-.' |  | |(   ) | (.-' 
   `--'`-'`-`-'  `-`--`-`-`-'`-`-'`--'

       doi:10.5281/zenodo.2634304
             Version {}
""".format(
            __version__
        )
    )


# For backwards-compatibility
say_hello = hello
