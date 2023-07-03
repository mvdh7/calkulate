import numpy as np
from scipy.optimize import least_squares
from matplotlib import pyplot as plt

# Density for H2SO4 at 298.15 K calculated with E-AIM Model Ic:
# http://www.aim.env.uea.ac.uk/aim/model1/model1c.php
molality = np.array(
    [
        0.05,
        0.1,
        0.2,
        0.5,
        1,
        1.5,
        2,
        2.5,
        3,
    ]
)
density = np.array(
    [
        1.00043,
        1.00361,
        1.00983,
        1.02785,
        1.05653,
        1.08378,
        1.10964,
        1.13410,
        1.15736,
    ]
)


# Generate fit
def density_fitter(coeffs, molality):
    return coeffs[0] + coeffs[1] * molality + coeffs[2] * molality**2


def lsq_density_fitter(coeffs, molality, density):
    return density_fitter(coeffs, molality) - density


density_fit = least_squares(
    lsq_density_fitter, [1.0, 0.05, 0.0], args=(molality, density)
)


fig, ax = plt.subplots(dpi=300)
ax.scatter(molality, density)
fx = np.linspace(np.min(molality), np.max(molality), num=500)
ax.plot(fx, density_fitter(density_fit["x"], fx))
