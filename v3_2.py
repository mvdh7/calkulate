from calkulate import v3_2 as calk

# Define metadata
metadata = calk.prepare_metadata(
    {
        "file_name": "tests/data/seawater-CRM-144.dat",
        "salinity": 33.3,
        "analyte_volume": 0.1,
    }
)


# Long-winded
tt = calk.read_dat(metadata["file_name"])
tt = calk.prepare_titration(metadata, tt)
gran_1 = calk.gran_estimator(tt["mixture_mass"], tt["emf"], tt["temperature"])

# Wrapped
gran_2 = calk.get_gran_estimator(metadata)
