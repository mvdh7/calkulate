# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019  Matthew Paul Humphreys  (GNU GPLv3)

from numpy import float_

# Universal physical constants
R     = float_(  8.3144621) # Gas constant     / J /K /mol
F     = float_( 96.4853365) # Faraday constant / kC/mol
Tzero = float_(273.15     ) # 0 degC           / K

# Salinity to chlorinity ratio [WLD69]
S_Cl = float_(1.80655)

# Relative molecular masses / g/mol
RMM_B = float_(10.811) # Boron
RMM_F = float_(18.998) # Fluoride
