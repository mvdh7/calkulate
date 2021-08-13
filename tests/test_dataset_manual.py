import copy
import calkulate as calk, numpy as np, pandas as pd

# Import file
file_name = "tests/data/seawater-CRM-144.dat"
titrant_volume, emf, temperature = calk.io.read_dat(file_name)
titrant_mass = titrant_volume * calk.density.HCl_NaCl_25C_DSC07() * 1e-3
analyte_mass = 0.1  # kg
dic = 2121  # micromol/kg-solution
salinity = 34.1  # practical
nutrients = {
    "total_silicate": 100,
    "total_phosphate": 20,
    "total_ammonia": 12,
    "total_sulfide": 13,
    "total_alpha": 25,
    "total_beta": 25,
}  # all micromol/kg-solution
k_alpha = 1e-5
k_beta = 1e-6
alkalinity_certified = 2345  # micromol/kg-solution


# Compare with calk.titration functions
ctf_kwargs = dict(
    analyte_mass=analyte_mass,
    dic=dic,
    **nutrients,
    k_alpha=k_alpha,
    k_beta=k_beta,
)
titrant_molinity__ctf, analyte_mass__ctf = calk.titration.calibrate(
    file_name, salinity, alkalinity_certified, **ctf_kwargs
)
(
    alkalinity__ctf,
    emf0__ctf,
    pH__ctf,
    temperature0__ctf,
    analyte_mass__ctf,
    opt_result__ctf,
) = calk.titration.solve(file_name, salinity, titrant_molinity__ctf, **ctf_kwargs)


# Test dataset tools
ds = pd.DataFrame({k: [v] for k, v in nutrients.items()})
ds["file_name"] = file_name
ds["analyte_mass"] = analyte_mass
ds["dic"] = dic
ds["salinity"] = salinity
ds["k_alpha"] = k_alpha
ds["k_beta"] = k_beta
ds["alkalinity_certified"] = alkalinity_certified
ds["file_good"] = True
ds2 = copy.deepcopy(ds)
ds3 = copy.deepcopy(ds)
ds4 = copy.deepcopy(ds)
ds5 = copy.deepcopy(ds)
ds5.loc[0, "file_name"] = "tests/data/PC_LIMS_Report-BATCH138-20200317-135120.txt"

# Calibrate a single row
ds_row = ds.loc[0]
titrant_molinity__ds_row = calk.dataset.calibrate_row(ds_row).to_numpy()[0]

# Calibrate the whole dataset
calk.dataset.calibrate(ds)
batches = calk.dataset.get_batches(ds)

# Solve a single row
ds_row = ds.loc[0]
solved__ds_row = calk.dataset.solve_row(ds_row)

# Solve the whole dataset
calk.dataset.solve(ds)
calk.dataset.calkulate(ds2)

# Use methods
ds3 = calk.dataset.Dataset(ds3)
ds3.calibrate()
ds3.solve()
ds4 = calk.dataset.Dataset(ds4)
ds4.calkulate()
ds5 = calk.Dataset(ds5)
ds5["read_dat_method"] = "pclims"
ds5.calkulate()


def test_calibrate_ds_row():
    """Do the calk.dataset functions give consistent calibration results?"""
    assert np.isclose(
        titrant_molinity__ds_row, titrant_molinity__ctf, rtol=0, atol=1e-12
    )
    assert np.isclose(
        titrant_molinity__ds_row, ds.loc[0, "titrant_molinity"], rtol=0, atol=1e-12
    )


def test_solve_ds_row():
    """Do the calk.dataset functions give consistent solver results?"""
    assert np.isclose(alkalinity__ctf, alkalinity_certified, rtol=0, atol=1e-6)
    assert np.isclose(alkalinity__ctf, solved__ds_row.alkalinity, rtol=0, atol=1e-12)
    assert np.isclose(emf0__ctf, solved__ds_row.emf0, rtol=0, atol=1e-12)
    assert np.isclose(pH__ctf, solved__ds_row.pH_initial, rtol=0, atol=1e-12)
    assert np.isclose(alkalinity__ctf, ds.loc[0, "alkalinity"], rtol=0, atol=1e-12)
    assert np.isclose(emf0__ctf, ds.loc[0, "emf0"], rtol=0, atol=1e-12)
    assert np.isclose(pH__ctf, ds.loc[0, "pH_initial"], rtol=0, atol=1e-12)


def test_methods():
    """Do the Dataset methods give consistent results?"""
    assert np.isclose(
        ds.loc[0, "alkalinity"], ds2.loc[0, "alkalinity"], rtol=0, atol=1e-6
    )
    assert np.isclose(
        ds2.loc[0, "alkalinity"], ds3.loc[0, "alkalinity"], rtol=0, atol=1e-6
    )
    assert np.isclose(
        ds3.loc[0, "alkalinity"], ds4.loc[0, "alkalinity"], rtol=0, atol=1e-6
    )
    assert np.isclose(ds5.loc[0, "alkalinity"], alkalinity_certified, rtol=0, atol=1e-6)


# test_calibrate_ds_row()
# test_solve_ds_row()
# test_methods()
