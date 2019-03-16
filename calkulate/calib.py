# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019  Matthew Paul Humphreys  (GNU GPLv3)

"""Calibrate acid concentration using reference material measurements."""

from scipy.optimize import least_squares as olsq
from . import solve


def halfGran(Macid, EMF, Tk, Msamp, AT_cert, XT, KX):
    """Calibrate acid concentration using half-Gran approach."""
    return olsq(lambda Cacid: solve.halfGran(Macid, EMF, Tk, Msamp, Cacid,
        *XT, *KX, suppress_warnings=True)[0] - AT_cert, 0.1, method='lm')


def MPH(Macid, EMF, Tk, Msamp, AT_cert, XT, KX):
    """Calibrate acid concentration using full least-squares approach."""
    return olsq(lambda Cacid: \
        solve.MPH(Macid, EMF, Tk, Msamp, Cacid, *XT, KX)['x'][0] - AT_cert,
        0.1, method='lm')


def DAA03(Macid, EMF, Tk, Msamp, AT_cert, XT, KX):
    """Calibrate acid concentration using Dickson et al. (2003) approach."""
    return olsq(lambda Cacid: \
        solve.DAA03(Macid, EMF, Tk, Msamp, Cacid, *XT, KX)['x'][0] - AT_cert,
        0.1, method='lm')
