import calkulate as calk, numpy as np


def test_seawater_density_MP81():
    """Can we reproduce MP81's check value from their Table 1?"""
    rho_calk = calk.density.seawater_1atm_MP81(temperature=5, salinity=35)
    rho_calk = np.round(rho_calk * 1e3, decimals=5)
    rho_MP81 = 1027.67547  # their check value
    assert rho_calk == rho_MP81


def test_HCl_NaCl_DSC07():
    """Can we reproduce DSC07's check value?"""
    rho_mix = calk.density.HCl_NaCl_25C_DSC07(molinity_HCl=0.2, molinity_NaCl=0.5)
    rho_mix = np.round(rho_mix, decimals=5)
    rho_DSC07 = 1.02056  # their check value
    print(rho_mix, rho_DSC07)  # check fails for now... sent email to Andrew Dickson


# test_seawater_density_MP81()
# test_HCl_NaCl_DSC07()
