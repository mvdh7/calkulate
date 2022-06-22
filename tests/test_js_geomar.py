import calkulate as calk, pandas as pd, numpy as np

# Import dataset and use dataset functions to calibrate and solve
tf = calk.read_csv("tests/data/js-geomar/js-geomar.csv")
tf.calkulate()
tf = pd.DataFrame(tf)

# Use single-titration functions to calibrate and solve
alkalinity_certified = 2193.07
tt_kwargs = dict(
    analyte_mass=0.0498537,
    dic=2048.36,
    total_phosphate=0.42,
    total_silicate=4,
    read_dat_method="pclims",
)
tt_molinity = calk.titration.calibrate(
    "tests/data/js-geomar/PC_LIMS_Report-20220518-124748.txt",
    33.425,
    alkalinity_certified,
    **tt_kwargs,
)[0]
tt_alkalinity = calk.titration.solve(
    "tests/data/js-geomar/PC_LIMS_Report-20220518-124748.txt",
    33.425,
    tt_molinity,
    **tt_kwargs,
)[0]


def test_js_titration_calibrate_solve():
    assert isinstance(tt_molinity, float)
    assert isinstance(tt_alkalinity, float)
    assert np.isclose(tt_alkalinity, alkalinity_certified)


# test_js_titration_calibrate_solve()
