import numpy as np
from . import constants


def _check_input_type(input_var, container):
    """Input type validity checker."""
    if container is None:
        assert not isinstance(input_var, str)
        return input_var
    else:
        assert input_var in container
        return container[input_var]


def gran_estimator(
    mixture_mass="mixture_mass",
    emf="emf",
    temperature="temperature",
    data=None,
    use_points=None,
):
    """Simple Gran-plot estimator function following DAA03 eq. 10."""
    mixture_mass = _check_input_type(mixture_mass, data)
    emf = _check_input_type(emf, data)
    temperature = _check_input_type(temperature, data)
    if use_points is None:
        G = np.full(np.size(emf), True)
    else:
        G = use_points
    temperature_K = temperature + constants.absolute_zero
    return mixture_mass[G] * np.exp(
        emf[G] * constants.faraday / (constants.ideal_gas * temperature_K[G])
    )
