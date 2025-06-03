# %%
import os

from calkulate import core
from calkulate.convert import amount_units, keys_cau
from calkulate.core import (
    keys_calibrate_emf,
    keys_calibrate_pH,
    keys_totals_ks,
    totals_ks,
)
from calkulate.meta import _get_kwargs_for
from calkulate.read.titrations import keys_read_dat, read_dat


file_name = "tests/data/seawater-CRM-144.dat"
kwargs = dict(analyte_mass=0.1)
salinity = 35
alkalinity_certified = 2300
