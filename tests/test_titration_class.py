import calkulate as calk

tt = calk.Titration(
    file_name="seawater-CRM-144.dat",
    file_path="tests/data/",
    analyte_mass=0.1,
    dic=2000,
)

tt.calibrate(2300)

x = tt.titration
y = tt.__dict__

# tt.plot_emf()
# tt.plot_pH()
# tt.plot_gran_emf0()
# tt.plot_gran_alkalinity()
x.plot.scatter("titrant_mass", "alkalinity_estimate")
