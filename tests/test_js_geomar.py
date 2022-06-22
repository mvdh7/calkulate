import calkulate as calk, pandas as pd


tf = calk.read_csv("tests/data/js-geomar/js-geomar.csv")
tf.calkulate()
tf = pd.DataFrame(tf)


test = calk.titration.calibrate(
    'tests/data/js-geomar/PC_LIMS_Report-20220518-124748.txt',
    35,2300, analyte_mass=0.05,read_dat_method='pclims')
