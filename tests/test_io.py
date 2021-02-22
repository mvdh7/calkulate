import numpy as np
import calkulate as calk

file_names = [
    "tests/data/PC_LIMS_Report-CRM1-20201211-115353.txt",
    "tests/data/PC_LIMS_Report-BATCH138-20200317-135120.txt",
    "tests/data/PC_LIMS_Report-SEA2-20200317-130328.txt",
]


def test_pclims_io():
    """Can PC LIMS Report files be imported correctly?"""
    for file_name in file_names:
        titrant_amount, measurement, temperature = calk.read_dat(
            file_name, method="pclims"
        )
        assert isinstance(titrant_amount, np.ndarray)
        assert isinstance(titrant_amount[0], float)
        assert isinstance(measurement, np.ndarray)
        assert isinstance(measurement[0], float)
        assert isinstance(temperature, np.ndarray)
        assert isinstance(temperature[0], float)
        assert (
            np.shape(titrant_amount) == np.shape(measurement) == np.shape(temperature)
        )


# test_pclims_io()
