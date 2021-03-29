# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019--2021  Matthew P. Humphreys  (GNU GPLv3)
"""Work with titration data in a file."""

import numpy as np, pandas as pd
from . import convert, core, default, density, interface, io


def get_dat_data(
    file_name,
    molinity_HCl=default.molinity_HCl,
    molinity_NaCl=default.molinity_NaCl,
    molinity_H2SO4=default.molinity_H2SO4,
    temperature_override=None,
    titrant=default.titrant,
    titrant_amount_unit=default.titrant_amount_unit,
    titrant_density=None,
    read_dat_method=default.read_dat_method,
    read_dat_kwargs={},
):
    """Import a dat file and convert titrant units to mass in kg."""
    # Import titration data file
    titrant_amount, emf, temperature = io.read_dat(
        file_name, method=read_dat_method, **read_dat_kwargs
    )
    if not pd.isnull(temperature_override):
        temperature = np.full_like(titrant_amount, temperature_override)
    # Get titrant mass
    if titrant_amount_unit == "ml":
        if not pd.isnull(titrant_density):
            titrant_mass = titrant_amount * titrant_density * 1e-3
        else:
            if titrant == "H2SO4":
                titrant_mass = (
                    titrant_amount * density.H2SO4_25C_EAIM(molinity_H2SO4) * 1e-3
                )
            else:
                titrant_mass = (
                    titrant_amount
                    * density.HCl_NaCl_25C_DSC07(
                        molinity_HCl=molinity_HCl, molinity_NaCl=molinity_NaCl,
                    )
                    * 1e-3
                )
    elif titrant_amount_unit == "g":
        titrant_mass = titrant_amount * 1e-3
    elif titrant_amount_unit == "kg":
        titrant_mass = titrant_amount
    else:
        print("titrant_amount_unit not recognised.")
    return titrant_mass, emf, temperature


def get_totals_k_constants(
    titrant_mass,
    temperature,
    analyte_mass,
    salinity,
    dic=0,
    total_alpha=0,
    total_ammonia=0,
    total_beta=0,
    total_phosphate=0,
    total_silicate=0,
    total_sulfide=0,
    total_borate=None,
    total_fluoride=None,
    total_sulfate=None,
    k_alpha=None,
    k_beta=None,
    k_ammonia=None,
    k_bisulfate=None,
    k_borate=None,
    k_carbonic_1=None,
    k_carbonic_2=None,
    k_fluoride=None,
    k_phosphoric_1=None,
    k_phosphoric_2=None,
    k_phosphoric_3=None,
    k_silicate=None,
    k_sulfide=None,
    k_water=None,
    opt_k_bisulfate=default.opt_k_bisulfate,
    opt_k_carbonic=default.opt_k_carbonic,
    opt_k_fluoride=default.opt_k_fluoride,
    opt_total_borate=default.opt_total_borate,
):
    # Get totals from PyCO2SYS
    totals, totals_pyco2 = interface.get_totals(
        salinity,
        dic=dic,
        total_alpha=total_alpha,
        total_beta=total_beta,
        total_ammonia=total_ammonia,
        total_phosphate=total_phosphate,
        total_silicate=total_silicate,
        total_sulfide=total_sulfide,
        total_borate=total_borate,
        total_fluoride=total_fluoride,
        total_sulfate=total_sulfate,
        opt_k_carbonic=opt_k_carbonic,
        opt_total_borate=opt_total_borate,
    )
    # Dilute totals with titrant
    totals = convert.dilute_totals(totals, titrant_mass, analyte_mass)
    totals_pyco2 = convert.dilute_totals_pyco2(totals_pyco2, titrant_mass, analyte_mass)
    # Get k_constants from PyCO2SYS
    k_constants = interface.get_k_constants(
        totals_pyco2,
        temperature,
        k_alpha=k_alpha,
        k_ammonia=k_ammonia,
        k_beta=k_beta,
        k_bisulfate=k_bisulfate,
        k_borate=k_borate,
        k_carbonic_1=k_carbonic_1,
        k_carbonic_2=k_carbonic_2,
        k_fluoride=k_fluoride,
        k_phosphoric_1=k_phosphoric_1,
        k_phosphoric_2=k_phosphoric_2,
        k_phosphoric_3=k_phosphoric_3,
        k_silicate=k_silicate,
        k_sulfide=k_sulfide,
        k_water=k_water,
        opt_k_bisulfate=opt_k_bisulfate,
        opt_k_carbonic=opt_k_carbonic,
        opt_k_fluoride=opt_k_fluoride,
        opt_total_borate=opt_total_borate,
    )
    return totals, k_constants


def prepare(
    file_name,
    salinity,
    analyte_mass=None,  # kg
    analyte_volume=None,  # ml
    dic=0,
    total_alpha=0,
    total_beta=0,
    total_ammonia=0,
    total_phosphate=0,
    total_silicate=0,
    total_sulfide=0,
    total_borate=None,
    total_fluoride=None,
    total_sulfate=None,
    k_alpha=None,
    k_ammonia=None,
    k_beta=None,
    k_bisulfate=None,
    k_borate=None,
    k_carbonic_1=None,
    k_carbonic_2=None,
    k_fluoride=None,
    k_phosphoric_1=None,
    k_phosphoric_2=None,
    k_phosphoric_3=None,
    k_silicate=None,
    k_sulfide=None,
    k_water=None,
    molinity_HCl=default.molinity_HCl,
    molinity_NaCl=default.molinity_NaCl,
    molinity_H2SO4=default.molinity_H2SO4,
    temperature_override=None,
    titrant=default.titrant,
    titrant_amount_unit=default.titrant_amount_unit,
    titrant_density=None,
    opt_k_bisulfate=default.opt_k_bisulfate,
    opt_k_carbonic=default.opt_k_carbonic,
    opt_k_fluoride=default.opt_k_fluoride,
    opt_total_borate=default.opt_total_borate,
    read_dat_kwargs={},
    read_dat_method=default.read_dat_method,
):
    """Prepare a titration data file for calibration and/or solving."""
    titrant_mass, emf, temperature = get_dat_data(
        file_name,
        molinity_HCl=molinity_HCl,
        molinity_NaCl=molinity_NaCl,
        molinity_H2SO4=molinity_H2SO4,
        temperature_override=temperature_override,
        titrant=titrant,
        titrant_density=titrant_density,
        titrant_amount_unit=titrant_amount_unit,
        read_dat_method=read_dat_method,
        read_dat_kwargs=read_dat_kwargs,
    )
    if pd.isnull(analyte_mass):
        analyte_mass = (
            analyte_volume
            * density.seawater_1atm_MP81(temperature=temperature[0], salinity=salinity,)
            * 1e-3
        )
    totals, k_constants = get_totals_k_constants(
        titrant_mass,
        temperature,
        analyte_mass,
        salinity,
        dic=dic,
        total_alpha=total_alpha,
        total_ammonia=total_ammonia,
        total_beta=total_beta,
        total_phosphate=total_phosphate,
        total_silicate=total_silicate,
        total_sulfide=total_sulfide,
        total_borate=total_borate,
        total_fluoride=total_fluoride,
        total_sulfate=total_sulfate,
        k_alpha=k_alpha,
        k_beta=k_beta,
        k_ammonia=k_ammonia,
        k_bisulfate=k_bisulfate,
        k_borate=k_borate,
        k_carbonic_1=k_carbonic_1,
        k_carbonic_2=k_carbonic_2,
        k_fluoride=k_fluoride,
        k_phosphoric_1=k_phosphoric_1,
        k_phosphoric_2=k_phosphoric_2,
        k_phosphoric_3=k_phosphoric_3,
        k_silicate=k_silicate,
        k_sulfide=k_sulfide,
        k_water=k_water,
        opt_k_bisulfate=opt_k_bisulfate,
        opt_k_carbonic=opt_k_carbonic,
        opt_k_fluoride=opt_k_fluoride,
        opt_total_borate=opt_total_borate,
    )
    return titrant_mass, emf, temperature, analyte_mass, totals, k_constants


def calibrate(
    file_name,
    salinity,
    alkalinity_certified,
    analyte_total_sulfate=None,
    titrant=default.titrant,
    titrant_molinity_guess=None,
    pH_range=default.pH_range,
    least_squares_kwargs=default.least_squares_kwargs,
    **prepare_kwargs,
):
    """Calibrate titrant_molinity for a titration file given alkalinity_certified."""
    titrant_mass, emf, temperature, analyte_mass, totals, k_constants = prepare(
        file_name, salinity, titrant=titrant, **prepare_kwargs
    )
    solver_kwargs = {
        "pH_range": pH_range,
        "least_squares_kwargs": least_squares_kwargs,
    }
    if titrant == "H2SO4":
        assert analyte_total_sulfate is not None
        titrant_molinity = core.calibrate_H2SO4(
            alkalinity_certified,
            titrant_mass,
            emf,
            temperature,
            analyte_mass,
            analyte_total_sulfate * 1e-6,
            salinity,
            totals,
            k_constants,
            titrant_molinity_guess=titrant_molinity_guess,
            least_squares_kwargs=least_squares_kwargs,
            solver_kwargs=solver_kwargs,
        )["x"][0]
    else:
        titrant_molinity = core.calibrate(
            alkalinity_certified,
            titrant_mass,
            emf,
            temperature,
            analyte_mass,
            totals,
            k_constants,
            titrant_molinity_guess=titrant_molinity_guess,
            least_squares_kwargs=least_squares_kwargs,
            solver_kwargs=solver_kwargs,
        )["x"][0]
    return titrant_molinity, analyte_mass


def solve(
    file_name,
    salinity,
    titrant_molinity,
    analyte_total_sulfate=None,
    titrant=default.titrant,
    pH_range=default.pH_range,
    least_squares_kwargs=default.least_squares_kwargs,
    **prepare_kwargs,
):
    """Solve alkalinity, EMF0 and initial pH for a titration file given titrant_molinity.
    
    Results in micromol/kg-solution, mV, and on the Free scale.
    """
    titrant_mass, emf, temperature, analyte_mass, totals, k_constants = prepare(
        file_name, salinity, titrant=titrant, **prepare_kwargs
    )
    if titrant == "H2SO4":
        # Update sulfate dilution by titrant and its consequences
        assert analyte_total_sulfate is not None
        totals["total_sulfate"] = (
            1e-6 * analyte_total_sulfate * analyte_mass
            + titrant_molinity * titrant_mass
        ) / (analyte_mass + titrant_mass)
        totals_pyco2 = convert.totals_to_pyco2(totals, salinity)
        k_constants = interface.get_k_constants(totals_pyco2, temperature)
        # Solve for alkalinity and EMF0
        opt_result = core.solve_emf_complete_H2SO4(
            titrant_molinity,
            titrant_mass,
            emf,
            temperature,
            analyte_mass,
            totals,
            k_constants,
            least_squares_kwargs=least_squares_kwargs,
            pH_range=pH_range,
        )
    else:
        opt_result = core.solve_emf_complete(
            titrant_molinity,
            titrant_mass,
            emf,
            temperature,
            analyte_mass,
            totals,
            k_constants,
            least_squares_kwargs=least_squares_kwargs,
            pH_range=pH_range,
        )
    alkalinity, emf0 = opt_result["x"]
    # Calculate initial pH
    pH_initial = convert.emf_to_pH(emf[0], emf0, temperature[0])
    return alkalinity * 1e6, emf0, pH_initial, temperature[0], analyte_mass
