# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019--2021  Matthew P. Humphreys  (GNU GPLv3)
"""Work with titration data in a file."""

import numpy as np, pandas as pd
from scipy.stats import linregress
import PyCO2SYS as pyco2
from . import convert, core, default, density, interface, io, plot


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
                        molinity_HCl=molinity_HCl,
                        molinity_NaCl=molinity_NaCl,
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
    opt_pH_scale=default.opt_pH_scale,
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
        opt_pH_scale=opt_pH_scale,
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
    opt_pH_scale=default.opt_pH_scale,
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
            * density.seawater_1atm_MP81(
                temperature=temperature[0],
                salinity=salinity,
            )
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
        opt_pH_scale=opt_pH_scale,
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
    emf0_guess=None,
    **prepare_kwargs,
):
    """Calibrate titrant_molinity for a titration file given alkalinity_certified."""
    titrant_mass, emf, temperature, analyte_mass, totals, k_constants = prepare(
        file_name, salinity, titrant=titrant, **prepare_kwargs
    )
    solver_kwargs = {
        "pH_range": pH_range,
        "least_squares_kwargs": least_squares_kwargs,
        "emf0_guess": emf0_guess,
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
    emf0_guess=None,
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
            emf0_guess=emf0_guess,
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
            emf0_guess=emf0_guess,
        )
    alkalinity, emf0 = opt_result["x"]
    # Calculate initial pH
    pH_initial = convert.emf_to_pH(emf[0], emf0, temperature[0])
    return alkalinity * 1e6, emf0, pH_initial, temperature[0], analyte_mass, opt_result


class Titration:
    def __init__(
        self,
        file_name="",
        file_path="",
        salinity=35,
        analyte_mass=None,
        analyte_volume=None,
        **prepare_kwargs,
    ):
        assert analyte_mass is not None or analyte_volume is not None
        self.file_name = file_name
        self.file_path = file_path
        self.salinity = salinity
        self.prepare_kwargs = prepare_kwargs.copy()
        self.prepare_kwargs["analyte_volume"] = analyte_volume
        self.prepare_kwargs["analyte_mass"] = analyte_mass
        if analyte_volume is not None:
            self.analyte_volume = analyte_volume
        (
            titrant_mass,
            emf,
            temperature,
            self.analyte_mass,
            totals,
            k_constants,
        ) = prepare(
            file_path + file_name,
            salinity,
            analyte_mass=analyte_mass,
            analyte_volume=analyte_volume,
            **prepare_kwargs,
        )
        self.titration = pd.DataFrame(
            {
                "titrant_mass": titrant_mass,
                "emf": emf,
                "temperature": temperature,
            }
        )
        self.titration["dilution_factor"] = convert.get_dilution_factor(
            titrant_mass, self.analyte_mass
        )
        for k, v in totals.items():
            self.titration[k] = v
        for k, v in k_constants.items():
            self.titration[k] = v
        self.calibrated = False
        self.solved = False

    def calibrate(self, alkalinity_certified, **calibrate_kwargs):
        self.alkalinity_certified = alkalinity_certified
        self.calibrate_kwargs = calibrate_kwargs.copy()
        self.titrant_molinity = calibrate(
            self.file_path + self.file_name,
            self.salinity,
            alkalinity_certified,
            **self.calibrate_kwargs,
            **self.prepare_kwargs,
        )[0]
        self.gran_guesses()
        pH_range = (
            self.calibrate_kwargs["pH_range"]
            if "pH_range" in self.calibrate_kwargs
            else default.pH_range
        )
        self.titration["G_final"] = (self.titration.pH_gran >= pH_range[0]) & (
            self.titration.pH_gran <= pH_range[1]
        )
        self.calibrated = True
        self.solve()

    def set_titrant_molinity(self, titrant_molinity):
        self.titrant_molinity = titrant_molinity
        self.gran_guesses()
        self.calibrated = False

    @np.errstate(invalid="ignore")
    def gran_guesses(self, emf0_guess=None):
        # Get simple Gran-plot estimator
        st = self.titration
        st["mixture_mass"] = st.titrant_mass + self.analyte_mass
        st["gran_estimates"] = core.gran_estimator(
            st.mixture_mass, st.emf, st.temperature
        )
        # Select which data points to use for first guesses
        st["G_gran"] = (st.gran_estimates >= 0.1 * np.max(st.gran_estimates)) & (
            st.gran_estimates <= 0.9 * np.max(st.gran_estimates)
        )
        # Make first guesses
        (
            self.alkalinity_gran,
            self.gran_slope,
            self.gran_intercept,
        ) = core.gran_guess_alkalinity(
            st.titrant_mass[st.G_gran],
            st.gran_estimates[st.G_gran],
            self.analyte_mass,
            self.titrant_molinity,
        )
        titrant = (
            self.prepare_kwargs["titrant"]
            if "titrant" in self.prepare_kwargs
            else default.titrant
        )
        if titrant == "H2SO4":
            self.alkalinity_gran *= 2
        if emf0_guess is None:
            st["emf0_gran"] = core.gran_guesses_emf0(
                st.titrant_mass,
                st.emf,
                st.temperature,
                self.analyte_mass,
                self.titrant_molinity,
                alkalinity_guess=self.alkalinity_gran,
                titrant=titrant,
                HF=0,
                HSO4=0,
            )
            self.emf0_gran = np.mean(st.emf0_gran[st.G_gran])
        else:
            self.emf0_gran = emf0_guess
        st["pH_gran"] = convert.emf_to_pH(st.emf, self.emf0_gran, st.temperature)
        self.alkalinity_gran *= 1e6

    def solve(self, titrant_molinity=None, **solve_kwargs):
        if titrant_molinity is not None:
            self.set_titrant_molinity(titrant_molinity)
        self.solve_kwargs = solve_kwargs.copy()
        pH_range = (
            self.solve_kwargs["pH_range"]
            if "pH_range" in self.solve_kwargs
            else default.pH_range
        )
        self.titration["G_final"] = (self.titration.pH_gran >= pH_range[0]) & (
            self.titration.pH_gran <= pH_range[1]
        )
        (
            self.alkalinity,
            self.emf0,
            self.pH_initial,
            self.pH_initial_temperature,
            self.analyte_mass,
            self.opt_result,
        ) = solve(
            self.file_path + self.file_name,
            self.salinity,
            self.titrant_molinity,
            **self.solve_kwargs,
            **self.prepare_kwargs,
        )
        self.titration["pH"] = convert.emf_to_pH(
            self.titration.emf, self.emf0, self.titration.temperature
        )
        self.do_CO2SYS()
        self.solved = True

    def do_CO2SYS(self):
        st = self.titration
        totals = {
            k: st[k].to_numpy() * 1e6 if k in st else 0
            for k in [
                "total_sulfate",
                "total_fluoride",
                "total_borate",
                "total_phosphate",
                "total_silicate",
                "total_ammonia",
                "total_sulfide",
                "total_alpha",
                "total_beta",
            ]
        }
        k_constants = {k: st[k].to_numpy() for k in st.columns if k.startswith("k_")}
        results = pyco2.sys(
            par1=st.dic.to_numpy() * 1e6,
            par2=st.pH.to_numpy(),
            par1_type=2,
            par2_type=3,
            opt_pH_scale=3,
            salinity=self.salinity,
            temperature=st.temperature,
            **totals,
            **k_constants,
        )
        for co2sysvar in [
            "alkalinity",
            "HCO3",
            "CO3",
            "BOH4",
            "PO4",
            "HPO4",
            "H3PO4",
            "HSO4",
            "HF",
            "H3SiO4",
            "OH",
        ]:
            st[co2sysvar] = results[co2sysvar] * 1e-6
        st["H"] = 10 ** -st.pH
        st["alk_alpha"] = results["alkalinity_alpha"] * 1e-6
        st["alk_beta"] = results["alkalinity_beta"] * 1e-6
        st["alkalinity_estimate"] = (
            st.alkalinity
            + st.titrant_mass
            * self.titrant_molinity
            / (st.titrant_mass + self.analyte_mass)
        ) / st.dilution_factor

    def calkulate(self, alkalinity_certified, calibrate_kwargs={}, solve_kwargs={}):
        self.calibrate(alkalinity_certified, **calibrate_kwargs)
        self.solve(**solve_kwargs)

    # Plotting functions
    plot_emf = plot.titration.emf
    plot_pH = plot.titration.pH
    plot_gran_emf0 = plot.titration.gran_emf0
    plot_gran_alkalinity = plot.titration.gran_alkalinity
    plot_alkalinity = plot.titration.alkalinity
    plot_components = plot.titration.components
