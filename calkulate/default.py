# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019--2021  Matthew P. Humphreys  (GNU GPLv3)
"""Set default values."""


dic = 0  # micromol / kg-solution
dpi = 300  # resolution of figures
opt_gas_constant = 3  # for PyCO2SYS
opt_k_bisulfate = 1  # for PyCO2SYS
opt_k_carbonic = 16  # for PyCO2SYS
opt_k_fluoride = 1  # for PyCO2SYS
opt_pH_scale = 3  # for PyCO2SYS
opt_total_borate = 1  # for PyCO2SYS
pressure = 0  # in-water pressure in dbar
pH_range = (3, 4)  # pH range to use for solving alkalinity
least_squares_kwargs = dict(method="lm", gtol=1e-12, xtol=1e-12)
molinity_HCl = 0.1  # for HCl acid density
molinity_NaCl = 0.6  # for HCl acid density
molinity_H2SO4 = 0.1  # for H2SO4 acid density
read_dat_method = "genfromtxt"
salinity = 35
total_ammonia = 0
total_phosphate = 0
total_silicate = 0
total_sulfide = 0
titrant = "HCl"
titrant_amount_unit = "ml"
titrant_molinity_guess = 0.1  # mol / kg-solution
verbose = False
zlp = 4.5  # pK of 'zero level of protons' [WZK07]
