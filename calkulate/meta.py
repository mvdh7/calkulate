# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019--2025  Matthew P. Humphreys  (GNU GPLv3)
"""Package metadata and metafunctions."""

import inspect

import pandas as pd


_authorlist = ["Humphreys, Matthew P.", "Matthews, Ruth S."]
__author__ = " and ".join(_authorlist)
__version__ = "23.7.0"
__year__ = 2025


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
""".format(__version__, __year__)
    )


def _get_kwarg_keys(func):
    params = inspect.signature(func).parameters
    return {p for p in params if params[p].default != inspect.Parameter.empty}


def _get_kwargs_for(keys, kwargs, row=None):
    # Start by getting any kwargs that are in the keys
    kwargs_for = {k: v for k, v in kwargs.items() if k in keys}
    # Overwrite / append any that (also) have values in the row
    if row is not None:
        for k in keys:
            if k in row and pd.notnull(row[k]):
                kwargs_for[k] = row[k]
    return kwargs_for
