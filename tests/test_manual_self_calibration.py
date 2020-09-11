import calkulate as calk, PyCO2SYS as pyco2, numpy as np, pandas as pd

# Import data
tf = calk.read_csv("tests/data/titration_table.csv").calkulate()
x = tf.titration[0]

# tf.apply(_get_totals, axis=1)
# print(tf.titration[0])
# print(tf.titration[0].columns)

# # # Manually prepare inputs for low-level functions
# # titrant = {
# #     "mass": tt.dat_dict["titrant_amount"] * 1e-3 / calk.density.HCl_NaCl_25C_DSC07(),
# # }
# # analyte = {
# #     "alkalinity_certified": tt.alkalinity_certified,
# #     "dic": tt.dic,
# #     "mass": tt.analyte_volume
# #     * 1e-3
# #     / calk.density.seawater_1atm_MP81(
# #         temperature=tt.dat_dict["mixture_temperature"][0], salinity=tt.salinity
# #     ),
# # }

# # # Evaluate PyCO2SYS intermediates
# # total_ammonia = 5
# # total_sulfide = 3
# # opt_k_carbonic = 16
# # totals = pyco2.salts.assemble(
# #     tt["salinity"],
# #     tt["total_silicate"],
# #     tt["total_phosphate"],
# #     total_ammonia,
# #     total_sulfide,
# #     opt_k_carbonic,
# #     tt["opt_total_borate"],
# # )
# # dilution_factor = calk.convert.dilution_factor(
# #     analyte["mass"], analyte["mass"] + titrant["mass"]
# # )
# # totals = {k: v * dilution_factor if k != "Sal" else v for k, v in totals.items()}
# # pressure = 0
# # opt_pH_scale = 3
# # opt_k_bisulfate = 1
# # opt_k_fluoride = 1
# # opt_gas_constant = 3
# # k_constants = pyco2.equilibria.assemble(
# #     tt.dat_dict["mixture_temperature"],
# #     pressure,
# #     totals,
# #     opt_pH_scale,
# #     opt_k_carbonic,
# #     opt_k_bisulfate,
# #     opt_k_fluoride,
# #     opt_gas_constant,
# # )

# # # Calibrate the sample with itself
# # calibrated = calk.solve.calibrate(
# #     titrant,
# #     analyte,
# #     tt.dat_dict["mixture_measurement"],
# #     tt.dat_dict["mixture_temperature"],
# #     totals,
# #     k_constants,
# # )
# # titrant["molinity"] = calibrated["x"][0]

# # # Solve the sample from its own calibration
# # solved = calk.solve.complete_emf(
# #     titrant,
# #     analyte,
# #     tt.dat_dict["mixture_measurement"],
# #     tt.dat_dict["mixture_temperature"],
# #     totals,
# #     k_constants,
# # )
# # alkalinity_calibrated = solved["x"][0] * 1e6


# # def test_manual_self_calibration():
# #     """Can we successfully calibrate a sample with itself?"""
# #     assert np.isclose(
# #         alkalinity_calibrated, analyte["alkalinity_certified"], rtol=0, atol=1e-8
# #     )


# # test_manual_self_calibration()
