# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019  Matthew Paul Humphreys  (GNU GPLv3)

from scipy.optimize import least_squares as olsq
from . import solve

def halfGran(Macid, emf, tempK, Msamp, AT_cert, XT, KXF):
    """Calibrate acid concentration using half-Gran method."""
    return olsq(lambda Cacid: solve.halfGran(Macid, emf, tempK, Msamp, Cacid,
        XT, KXF, suppress_warnings=True)[0] - AT_cert, 0.1, method='lm')

def complete(Macid, emf, tempK, Msamp, AT_cert, XT, KXF):
    """Calibrate acid concentration using Complete Calculation."""
    return olsq(lambda Cacid: solve.complete(Macid, emf, tempK, Msamp, Cacid, 
        XT, KXF)['x'][0] - AT_cert, 0.1, method='lm')

def DAA03(Macid, emf, tempK, Msamp, AT_cert, XT, KXF):
    """Calibrate acid concentration following Dickson et al. (2003)."""
    return olsq(lambda Cacid: solve.DAA03(Macid, emf, tempK, Msamp, Cacid, 
        XT, KXF)['x'][0] - AT_cert, 0.1, method='lm')
