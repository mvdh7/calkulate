from calkulate import v3_2 as calk

tt = calk.read_dat("tests/data/seawater-CRM-144.dat")

test = calk.gran_estimator(
    mixture_mass="titrant_amount",
    emf="measurement",
    # temperature="temperature",
    data=tt,
    # use_points=None,
)

test2 = calk.gran_estimator(
    mixture_mass=tt.titrant_amount,
    emf=tt.measurement,
    temperature=tt.temperature,
    # data=None,
    # use_points=None,
)
