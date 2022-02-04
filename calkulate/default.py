# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019--2022  Matthew P. Humphreys  (GNU GPLv3)
"""Set default values."""


dic = 0  # micromol / kg-solution
dpi = 300  # resolution of figures
fCO2_air = 450  # for DIC loss modelling, µatm
least_squares_kwargs = dict(method="lm", gtol=1e-12, xtol=1e-12)
molinity_H2SO4 = 0.1  # for H2SO4 acid density
molinity_HCl = 0.1  # for HCl acid density
molinity_NaCl = 0.6  # for HCl acid density
opt_gas_constant = 3  # for PyCO2SYS
opt_k_bisulfate = 1  # for PyCO2SYS
opt_k_carbonic = 16  # for PyCO2SYS
opt_k_fluoride = 1  # for PyCO2SYS
opt_pH_scale = 3  # for PyCO2SYS
opt_total_borate = 1  # for PyCO2SYS
pH_range = (3, 4)  # pH range to use for solving alkalinity
pressure = 0  # in-water pressure in dbar
read_dat_method = "genfromtxt"
salinity = 35  # practical
split_pH = 5.5  # for DIC loss modelling
titrant = "HCl"
titrant_amount_unit = "ml"
titrant_molinity_guess = 0.1  # mol / kg-solution
total_ammonia = 0  # µmol/kg-sw
total_phosphate = 0  # µmol/kg-sw
total_silicate = 0  # µmol/kg-sw
total_sulfide = 0  # µmol/kg-sw
verbose = False
zlp = 4.5  # pK of 'zero level of protons' [WZK07]
