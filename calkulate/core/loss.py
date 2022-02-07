# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019--2022  Matthew P. Humphreys  (GNU GPLv3)
"""Model DIC loss during titrations."""

import numpy as np
from scipy import interpolate, optimize
from .. import convert, default


def dic_loss_model_fitted(
    k_dic_loss, i_titrant_mass, i_delta_fCO2, dic_start, analyte_mass, step_tm
):
    """Calculate DIC loss for the fitted section of the titration."""
    i_dic = np.full_like(i_delta_fCO2, np.nan)
    i_dic[0] = dic_start * 1e6
    for i in range(len(i_dic) - 1):
        dilute = convert.get_dilution_factor(step_tm, analyte_mass + i_titrant_mass[i])
        i_dic[i + 1] = i_dic[i] * dilute - k_dic_loss * i_delta_fCO2[i] * step_tm
    return i_dic


def _lsqfun_dic_loss_model(
    k_dic_loss,
    i_titrant_mass,
    i_delta_fCO2,
    dic_start,
    analyte_mass,
    step_tm,
    i_dic_loss,
):
    return (
        dic_loss_model_fitted(
            k_dic_loss, i_titrant_mass, i_delta_fCO2, dic_start, analyte_mass, step_tm
        )
        - i_dic_loss
    )


def get_fCO2_from_dic_pH(dic, pH, k0, k1, k2):
    """Calculate fCO2 from DIC and pH."""
    h = 10 ** -pH
    CO2aq = dic / (1 + k1 / h + k1 * k2 / h ** 2)
    fCO2 = CO2aq / k0
    return fCO2


def dic_loss_model_future(
    k_dic_loss,
    f_titrant_mass,
    f_pH,
    f_k0,
    f_k1,
    f_k2,
    delta_fCO2_start,
    dic_start,
    pH_start,
    k0_start,
    k1_start,
    k2_start,
    analyte_mass,
    step_tm,
    fCO2_air=default.fCO2_air,
):
    """Forecast future DIC loss, starting from the end of the fitted section."""
    # Prepare empty arrays for model
    f_delta_fCO2 = np.full_like(f_titrant_mass, np.nan)
    f_dic = np.full_like(f_titrant_mass, np.nan)
    # Get first values in the forecast arrays
    dilute = convert.get_dilution_factor(
        step_tm, analyte_mass + f_titrant_mass[0] - step_tm
    )
    f_dic[0] = dic_start * dilute - k_dic_loss * delta_fCO2_start * step_tm
    f_delta_fCO2[0] = (
        get_fCO2_from_dic_pH(f_dic[0], pH_start, k0_start, k1_start, k2_start)
        - fCO2_air
    )
    # Forecast future DIC loss
    for i in range(len(f_dic) - 1):
        dilute = convert.get_dilution_factor(step_tm, analyte_mass + f_titrant_mass[i])
        f_dic[i + 1] = f_dic[i] * dilute - k_dic_loss * f_delta_fCO2[i] * step_tm
        f_delta_fCO2[i + 1] = (
            get_fCO2_from_dic_pH(
                f_dic[i + 1], f_pH[i + 1], f_k0[i + 1], f_k1[i + 1], f_k2[i + 1]
            )
            - fCO2_air
        )
    return f_dic, f_delta_fCO2


def get_dic_loss_hires(
    titrant_mass,
    pH,
    dic_loss,
    fCO2_loss,
    k_CO2,
    k_carbonic_1,
    k_carbonic_2,
    analyte_mass,
    dic_start,
    fCO2_air=default.fCO2_air,
    split_pH=default.split_pH,
):
    """Fit and forecast high-resolution DIC loss model."""
    # Get delta-fCO2
    delta_fCO2_loss = fCO2_loss - fCO2_air
    # Use titrant_mass as proxy for titration time: generate high-resolution arrays
    step_tm = np.median(np.diff(titrant_mass)) / 10
    a_titrant_mass = np.arange(0, np.max(titrant_mass), step_tm)
    a_pH = interpolate.pchip_interpolate(titrant_mass, pH, a_titrant_mass)
    i_titrant_mass = a_titrant_mass[a_pH >= split_pH]
    f_titrant_mass = a_titrant_mass[a_pH < split_pH]
    # Interpolate other properties to the high-resolution titrant_mass
    i_dic_loss = interpolate.pchip_interpolate(titrant_mass, dic_loss, i_titrant_mass)
    i_delta_fCO2 = interpolate.pchip_interpolate(
        titrant_mass, delta_fCO2_loss, i_titrant_mass
    )
    f_pH = interpolate.pchip_interpolate(titrant_mass, pH, f_titrant_mass)
    f_k0 = interpolate.pchip_interpolate(titrant_mass, k_CO2, f_titrant_mass)
    f_k1 = interpolate.pchip_interpolate(titrant_mass, k_carbonic_1, f_titrant_mass)
    f_k2 = interpolate.pchip_interpolate(titrant_mass, k_carbonic_2, f_titrant_mass)
    # Get mid-way start-points for forecasting
    pH_start = interpolate.pchip_interpolate(titrant_mass, pH, i_titrant_mass[-1])
    k0_start = interpolate.pchip_interpolate(titrant_mass, k_CO2, i_titrant_mass[-1])
    k1_start = interpolate.pchip_interpolate(
        titrant_mass, k_carbonic_1, i_titrant_mass[-1]
    )
    k2_start = interpolate.pchip_interpolate(
        titrant_mass, k_carbonic_2, i_titrant_mass[-1]
    )
    # Find best-fit k_dic_loss
    k_dic_loss_opt_result = optimize.least_squares(
        _lsqfun_dic_loss_model,
        1.0,
        args=(
            i_titrant_mass,
            i_delta_fCO2,
            dic_start,
            analyte_mass,
            step_tm,
            i_dic_loss,
        ),
    )
    k_dic_loss = k_dic_loss_opt_result["x"]
    # Calculate DIC from k_dic_loss in the fitted region
    i_dic = dic_loss_model_fitted(
        k_dic_loss, i_titrant_mass, i_delta_fCO2, dic_start, analyte_mass, step_tm
    )
    # Forecast future DIC loss
    f_dic, f_delta_fCO2 = dic_loss_model_future(
        k_dic_loss,
        f_titrant_mass,
        f_pH,
        f_k0,
        f_k1,
        f_k2,
        i_delta_fCO2[-1],
        i_dic[-1],
        pH_start,
        k0_start,
        k1_start,
        k2_start,
        analyte_mass,
        step_tm,
        fCO2_air=fCO2_air,
    )
    return k_dic_loss[0], {
        "titrant_mass": a_titrant_mass,
        "pH": a_pH,
        "dic": np.concatenate((i_dic, f_dic)),
        "delta_fCO2": np.concatenate((i_delta_fCO2, f_delta_fCO2)),
    }


def get_dic_loss(
    titrant_mass,
    pH,
    dic_loss,
    fCO2_loss,
    k_CO2,
    k_carbonic_1,
    k_carbonic_2,
    analyte_mass,
    dic_start,
    fCO2_air=default.fCO2_air,
    split_pH=default.split_pH,
):
    """Get final DIC loss values at the titration points to go in the titration df."""
    k_dic_loss, loss_hires = get_dic_loss_hires(
        titrant_mass,
        pH,
        dic_loss,
        fCO2_loss,
        k_CO2,
        k_carbonic_1,
        k_carbonic_2,
        analyte_mass,
        dic_start,
        fCO2_air=fCO2_air,
        split_pH=split_pH,
    )
    dic_loss_modelled = (
        interpolate.pchip_interpolate(
            loss_hires["titrant_mass"],
            loss_hires["dic"],
            titrant_mass,
        )
        * 1e-6
    )
    fCO2_loss_modelled = (
        interpolate.pchip_interpolate(
            loss_hires["titrant_mass"],
            loss_hires["delta_fCO2"],
            titrant_mass,
        )
        + fCO2_air
    )
    dic_loss_fitted = pH >= split_pH
    return k_dic_loss, dic_loss_modelled, fCO2_loss_modelled, dic_loss_fitted
