import calkulate as calk, numpy as np

# Set up random number generator for repeatability
rng = np.random.default_rng(7)


def test_emf_h_conversion():
    """Can EMF and [H+] be accurately interconverted?"""
    tpts = 100
    emf0 = np.append(
        np.array([0.0, 100.0, 0.0]), rng.normal(size=tpts, scale=100, loc=0)
    )
    emf_original = np.append(
        np.array([0.0, 0.0, 100.0]), rng.normal(size=tpts, scale=200, loc=0)
    )
    temperature = np.append(
        np.array([25.0, 25.0, 25.0]), rng.normal(size=tpts, scale=10, loc=25)
    )
    h = calk.convert.emf_to_h(emf_original, emf0, temperature=temperature)
    assert np.all(h > 0)
    emf_converted = calk.convert.h_to_emf(h, emf0, temperature=temperature)
    emf_diff = emf_converted - emf_original
    assert np.all(np.isclose(emf_diff, 0, rtol=0, atol=1e-12))


def test_emf_pH_conversion():
    """Can EMF and pH be accurately interconverted?"""
    tpts = 100
    emf0 = np.append(
        np.array([0.0, 100.0, 0.0]), rng.normal(size=tpts, scale=100, loc=0)
    )
    emf_original = np.append(
        np.array([0.0, 0.0, 100.0]), rng.normal(size=tpts, scale=200, loc=0)
    )
    temperature = np.append(
        np.array([25.0, 25.0, 25.0]), rng.normal(size=tpts, scale=10, loc=25)
    )
    pH = calk.convert.emf_to_pH(emf_original, emf0, temperature=temperature)
    emf_converted = calk.convert.pH_to_emf(pH, emf0, temperature=temperature)
    emf_diff = emf_converted - emf_original
    assert np.all(np.isclose(emf_diff, 0, rtol=0, atol=1e-12))


# test_emf_h_conversion()
# test_emf_pH_conversion()
