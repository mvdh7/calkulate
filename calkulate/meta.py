# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019--2021  Matthew P. Humphreys  (GNU GPLv3)
"""Package metadata."""

_authorlist = ["Humphreys, Matthew P.", "Matthews, Ruth S."]
__author__ = " and ".join(_authorlist)
__version__ = "23.2.1"


def hello():
    """Report the version number and DOI."""
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
