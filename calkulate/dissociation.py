# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew Paul Humphreys  (GNU GPLv3)
"""Evaluate stoichiometric equilibrium constants."""
from PyCO2SYS import assemble

def eqConstants(tempK, pSal, concTotals, WhichKs=10, WhoseKSO4=1, WhoseKF=1):
    """Assemble a dict of dissociation constants on the Free pH scale.
    Use PyCO2SYS functions.
    """
    eqConstants = {} # for Free scale dissociation constants
    econsts = assemble.inputs({
        'TEMPIN': tempK - 273.15,
        'PRESIN': 0,
        'pHScale': 3,
        'WhichKs': WhichKs,
        'WhoseKSO4': WhoseKSO4,
        'WhoseKF': WhoseKF,
        'TP': concTotals['P'],
        'TSi': concTotals['Si'],
        'SAL': pSal,
        'TF': concTotals['F'],
        'TS': concTotals['S'],
        })[0].values()
    (_, eqConstants['C1'], eqConstants['C2'], eqConstants['w'],
            eqConstants['B'], eqConstants['F'], eqConstants['S'],
            eqConstants['P1'], eqConstants['P2'], eqConstants['P3'],
            eqConstants['Si'], eqConstants['NH3'], eqConstants['H2S'], _) = \
        assemble.equilibria(*econsts)
    return eqConstants
