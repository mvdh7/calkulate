import numpy as np
import calkulate as calk


def test_emf_h_conversion():
    """Check that EMF and [H+] can be accurately interconverted."""
    tpts = 100
    emf0 = np.append(
        np.array([0.0, 100.0, 0.0]), np.random.normal(size=tpts, scale=100, loc=0)
    )
    emf_original = np.append(
        np.array([0.0, 0.0, 100.0]), np.random.normal(size=tpts, scale=200, loc=0)
    )
    temperature = np.append(
        np.array([25.0, 25.0, 25.0]), np.random.normal(size=tpts, scale=10, loc=25)
    )
    h = calk.convert.emf_to_h(emf_original, emf0, temperature)
    assert np.all(h > 0)
    emf_converted = calk.convert.h_to_emf(h, emf0, temperature)
    emf_diff = emf_converted - emf_original
    assert np.max(np.abs(emf_diff)) < 1e-12


test_emf_h_conversion()
