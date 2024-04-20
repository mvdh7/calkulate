# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019--2024  Matthew P. Humphreys  (GNU GPLv3)
"""Package metadata."""

_authorlist = ["Humphreys, Matthew P.", "Matthews, Ruth S."]
__author__ = " and ".join(_authorlist)
__version__ = "23.6.1"
__year__ = 2024


def hello():
    """Report the version number and DOI."""
    print(
        r"""
   .--.     . .         .      .      
  :         | |         |     _|_     
  |    .-.  | |.-. .  . | .-.  |  .-. 
  :   (   ) | |-.' |  | |(   ) | (.-' 
   `--'`-'`-`-'  `-`--`-`-`-'`-`-'`--'

    Humphreys, M.P. & Matthews, R.S.
       doi:10.5281/zenodo.2634304
           Version {} ({})
""".format(
            __version__, __year__
        )
    )
