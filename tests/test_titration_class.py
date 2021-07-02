import calkulate as calk

tt = calk.Titration(
    file_name="seawater-CRM-144.dat",
    file_path="tests/data/",
    analyte_mass=0.1,
    dic=3000,
)

tt.calibrate(2300)
tt.solve()
# tt.calkulate(2300)
