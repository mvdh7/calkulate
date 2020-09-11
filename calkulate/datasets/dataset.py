import pandas as pd


class Dataset(pd.DataFrame):
    from .get import (
        get_measurement_type,
        get_titrations,
        get_analyte_temperature,
        get_analyte_mass,
        get_titrant_density,
        get_titrant_mass,
        get_analyte_totals,
        get_titration_totals,
        get_totals,
        get_k_constants,
        calkulate,
    )
    from .quantify import solve
