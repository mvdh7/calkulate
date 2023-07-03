# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019--2023  Matthew P. Humphreys  (GNU GPLv3)
"""Work with titration data in a file."""

import numpy as np, pandas as pd
from scipy.stats import linregress
import PyCO2SYS as pyco2
from . import convert, core, default, density, interface, read, plot, simulate


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
    titrant_amount, emf, temperature = read.read_dat(
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
    dic=default.dic,
    total_alpha=0,
    total_ammonia=default.total_ammonia,
    total_beta=0,
    total_phosphate=default.total_phosphate,
    total_silicate=default.total_silicate,
    total_sulfide=default.total_sulfide,
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
    dic=default.dic,
    total_alpha=0,
    total_beta=0,
    total_ammonia=default.total_ammonia,
    total_phosphate=default.total_phosphate,
    total_silicate=default.total_silicate,
    total_sulfide=default.total_sulfide,
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
    """Calibrate ``titrant_molinity`` for a titration file given
    ``alkalinity_certified``.
    """
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
    """Solve alkalinity, EMF0 and initial pH for a titration file given the
    ``titrant_molinity``.

    Results in micromol/kg-solution, mV, and on the Free pH scale.
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
    """
    ``calk.Titration``
    ==================
    A structure containing all the information about one titration with methods for
    calibrating and solving and then visualising the results.

    Creating a ``Titration``
    ------------------------
    There are two ways to create a ``Titration``:

    (1) From a row in ``calk.Dataset`` using the ``to_Titration`` method (or
    equivalently from a ``pandas.DataFrame`` using the ``calk.dataset.to_Titration``
    function):

      >>> tt = ds.to_Titration(index)
      >>> tt = calk.dataset.to_Titration(ds, index)

    (2) By initalising a ``calk.Titration`` directly with a ``file_name`` for a
    titration data file.

      >>> tt = calk.Titration(file_name="path/to/file.dat", **kwargs)

    Alkalinity solving methods
    --------------------------
    ``calibrate``
        Find the best-fit ``titrant_molinity`` for a given ``alkalinity_certified``.
    ``set_titrant_molinity``
        Set the ``titrant_molinity`` if it is known independently.
    ``solve``
        Solve for ``alkalinity`` when ``titrant_molinity`` is known.
    ``calkulate``
        Run the ``calibrate`` and ``solve`` steps sequentially.

    Data visualisation methods
    --------------------------
    ``plot_emf``
        EMF throughout the titration.
    ``plot_pH``
        pH throughout the titration.
    ``plot_gran_emf0``
        The Gran-plot initial estimate of EMF0.
    ``plot_gran_alkalinity``
        The Gran-plot initial estimate of alkalinity.
    ``plot_alkalinity``
        The final alkalinity estimates across the titration.
    ``plot_components``
        Show how all equilibrating components of the solution change throughout the
        titration.

    Attributes
    ----------
    ``alkalinity`` : ``float``
        The solved alkalinity in micromol/kg-seawater.
    ``alkalinity_certified`` : ``float``
        The certified alkalinity value against which the sample was calibrated,
        if applicable, in micromol/kg-seawater.
    ``alkalinity_gran`` : ``float``
        The Gran-plot initial estimate of alkalinity in micromol/kg-seawater.
    ``analyte_mass`` : ``float``
        The mass of the analyte in kg.
    ``analyte_volume`` : ``float``
        The volume of the analyte in ml, if provided.
    ``calibrated`` : ``bool``
        Whether or not the ``Titration`` has been calibrated.
    ``dic`` : ``float``
        The dissolved inorganic carbon in micromol/kg-seawater.
    ``emf0`` : ``float``
        The solved EMF0 in mV.
    ``emf0_gran`` : ``float``
        The Gran-plot initial estimate of EMF0 in mV.
    ``file_name`` : ``str``
        The name of the titration data file.
    ``file_path`` : ``str``
        The path to the titration data file.
    ``file_prepare_kwargs`` : ``dict``
        The kwargs passed to prepare the titration data file for solving.
    ``opt_result`` : ``scipy.optimize._optimize.OptimizeResult``
        The raw output from the ``scipy.optimize.least_squares`` solver.
    ``pH_initial`` : ``float``
        The pH at the initial titration point on the free pH scale.
    ``pH_initial_temperature`` : ``float``
        The temperature at the initial titration point.
    ``pH_range`` : ``tuple`` or ``list``
        The (minimum, maximum) pH to use for solving the titration.
    ``salinity`` : ``float``
        The practical salinity of the analyte.
    ``solved`` : ``bool``
        Whether or not the titration has been solved.
    ``titrant`` : ``str``
        The type of titrant used in the titration.
    ``titrant_molinity`` : ``float``
        The molinity of the titrant in mol/kg-titrant.
    ``titration`` : ``pd.DataFrame``
        A table of everything that changes through the titration (e.g., amount of
        titrant, EMF, temperature, total salt concentrations and equilibrium constants).
    """

    def __init__(
        self,
        file_name=None,
        file_path="",
        salinity=None,
        analyte_mass=None,
        analyte_volume=None,
        file_prepare_kwargs={},
        simulate_alkalinity=None,
        simulate_kwargs={},
    ):
        assert (
            analyte_mass is not None or analyte_volume is not None
        ), "You must provide either analyte_mass [kg] or analyte_volume [ml]!"
        self.file_name = file_name
        self.file_path = file_path
        self.check_set_default("salinity", salinity)
        self.analyte_volume = analyte_volume
        if file_name is not None:
            self.file_prepare_kwargs = file_prepare_kwargs.copy()
            (
                titrant_mass,
                emf,
                temperature,
                self.analyte_mass,
                totals,
                k_constants,
            ) = prepare(
                self.file_path + self.file_name,
                self.salinity,
                analyte_mass=analyte_mass,
                analyte_volume=self.analyte_volume,
                **self.file_prepare_kwargs,
            )
            if "dic" in file_prepare_kwargs:
                self.dic = file_prepare_kwargs["dic"]
            else:
                self.dic = default.dic
        else:
            assert simulate_alkalinity is not None
            (
                titrant_mass,
                emf,
                temperature,
                self.analyte_mass,
                totals,
                k_constants,
            ) = simulate._titration(simulate_alkalinity, **simulate_kwargs)
            if "dic" in simulate_kwargs:
                self.dic = simulate_kwargs["dic"]
            else:
                self.dic = default.dic
        # Now do the processing that's independent of the data source
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
        # Put totals and k_constants into the titration df and get lists of their names
        self._totals = []
        self._k_constants = []
        for k, v in totals.items():
            self.titration[k] = v
            self._totals.append(k)
        for k, v in k_constants.items():
            self.titration[k] = v
            self._k_constants.append(k)
        self.calibrated = False
        self.solved = False

    def check_set_default(self, attr, value):
        """If value is None, check if attribute already exists.  If it does, do nothing,
        but if it doesn't, set it to the default.  If value is not None, set attribute
        to the new value, creating it if necessary.
        """
        if value is None:
            if not hasattr(self, attr):
                setattr(self, attr, getattr(default, attr))
        else:
            setattr(self, attr, value)

    def get_totals(self):
        return {k: self.titration[k].to_numpy() for k in self._totals}

    def get_k_constants(self):
        return {k: self.titration[k].to_numpy() for k in self._k_constants}

    def calibrate(
        self,
        alkalinity_certified,
        analyte_total_sulfate=None,
        emf0_guess=None,
        least_squares_kwargs=None,
        pH_range=None,
        titrant_molinity_guess=None,
        titrant=None,
    ):
        self.alkalinity_certified = alkalinity_certified
        self.analyte_total_sulfate = analyte_total_sulfate
        self.check_set_default("least_squares_kwargs", least_squares_kwargs)
        self.check_set_default("pH_range", pH_range)
        self.check_set_default("titrant", titrant)
        # Solve for titrant_molinity
        st = self.titration
        solver_kwargs = {
            "emf0_guess": emf0_guess,
            "least_squares_kwargs": least_squares_kwargs,
            "pH_range": self.pH_range,
        }
        if titrant == "H2SO4":
            assert self.analyte_total_sulfate is not None
            self.titrant_molinity = core.calibrate_H2SO4(
                self.alkalinity_certified,
                st.titrant_mass.to_numpy(),
                st.emf.to_numpy(),
                st.temperature.to_numpy(),
                self.analyte_mass,
                self.analyte_total_sulfate * 1e-6,
                self.salinity,
                self.get_totals(),
                self.get_k_constants(),
                titrant_molinity_guess=titrant_molinity_guess,
                least_squares_kwargs=least_squares_kwargs,
                solver_kwargs=solver_kwargs,
            )["x"][0]
        else:
            self.titrant_molinity = core.calibrate(
                self.alkalinity_certified,
                st.titrant_mass.to_numpy(),
                st.emf.to_numpy(),
                st.temperature.to_numpy(),
                self.analyte_mass,
                self.get_totals(),
                self.get_k_constants(),
                titrant_molinity_guess=titrant_molinity_guess,
                least_squares_kwargs=least_squares_kwargs,
                solver_kwargs=solver_kwargs,
            )["x"][0]
        # Get Gran-plot guesses with solved titrant molinity
        self.gran_guesses()
        self.calibrated = True

    def set_titrant_molinity(self, titrant_molinity):
        self.titrant_molinity = titrant_molinity
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
            alkalinity_gran,
            self.gran_slope,
            self.gran_intercept,
        ) = core.gran_guess_alkalinity(
            st.titrant_mass[st.G_gran],
            st.gran_estimates[st.G_gran],
            self.analyte_mass,
            self.titrant_molinity,
        )
        if self.titrant == "H2SO4":
            alkalinity_gran *= 2
        if emf0_guess is None:
            st["emf0_gran"] = core.gran_guesses_emf0(
                st.titrant_mass,
                st.emf,
                st.temperature,
                self.analyte_mass,
                self.titrant_molinity,
                alkalinity_guess=alkalinity_gran,
                titrant=self.titrant,
                HF=0,
                HSO4=0,
            )
            self.emf0_gran = np.mean(st.emf0_gran[st.G_gran])
        else:
            self.emf0_gran = emf0_guess
        st["pH_gran"] = convert.emf_to_pH(st.emf, self.emf0_gran, st.temperature)
        self.alkalinity_gran = alkalinity_gran * 1e6
        self.titration["G_final"] = (self.titration.pH_gran >= self.pH_range[0]) & (
            self.titration.pH_gran <= self.pH_range[1]
        )

    def solve(
        self,
        analyte_total_sulfate=None,
        emf0_guess=None,
        least_squares_kwargs=default.least_squares_kwargs,
        pH_range=None,
        titrant_molinity=None,
        titrant=None,
    ):
        if analyte_total_sulfate is not None:
            self.analyte_total_sulfate = analyte_total_sulfate
        if titrant_molinity is not None:
            self.set_titrant_molinity(titrant_molinity)
        self.check_set_default("pH_range", pH_range)
        self.check_set_default("titrant", titrant)
        self.gran_guesses()
        # Solve for total alkalinity and EMF0
        st = self.titration
        if self.titrant == "H2SO4":
            # Update sulfate dilution by titrant and its consequences
            assert self.analyte_total_sulfate is not None
            st["total_sulfate"] = (
                1e-6 * self.analyte_total_sulfate * self.analyte_mass
                + st.titrant_molinity * st.titrant_mass
            ) / (self.analyte_mass + st.titrant_mass)
            # Solve for alkalinity and EMF0
            self.opt_result = core.solve_emf_complete_H2SO4(
                self.titrant_molinity,
                st.titrant_mass.to_numpy(),
                st.emf.to_numpy(),
                st.temperature.to_numpy(),
                self.analyte_mass,
                self.get_totals(),
                self.get_k_constants(),
                least_squares_kwargs=least_squares_kwargs,
                pH_range=self.pH_range,
                emf0_guess=emf0_guess,
            )
        else:
            self.opt_result = core.solve_emf_complete(
                self.titrant_molinity,
                st.titrant_mass.to_numpy(),
                st.emf.to_numpy(),
                st.temperature.to_numpy(),
                self.analyte_mass,
                self.get_totals(),
                self.get_k_constants(),
                least_squares_kwargs=least_squares_kwargs,
                pH_range=self.pH_range,
                emf0_guess=emf0_guess,
            )
        self.alkalinity, self.emf0 = self.opt_result["x"]
        self.alkalinity *= 1e6
        # Calculate initial pH
        st["pH"] = convert.emf_to_pH(st.emf, self.emf0, st.temperature)
        self.pH_initial = st["pH"].iloc[0]
        self.pH_initial_temperature = st["temperature"].iloc[0]
        self.get_alkalinity_from_acid()
        self.do_CO2SYS()
        self.solved = True

    def get_alkalinity_from_acid(self):
        self.titration["alkalinity_from_acid"] = (
            1e-6 * self.alkalinity * self.analyte_mass
            - self.titration.titrant_mass * self.titrant_molinity
        ) / (self.analyte_mass + self.titration.titrant_mass)

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
            temperature=st.temperature.to_numpy(),
            **totals,
            **k_constants,
        )
        st["alkalinity_from_pH"] = results["alkalinity"] * 1e-6
        st["k_CO2"] = results["k_CO2"]
        for co2sysvar in [
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
        st["H"] = 10**-st.pH
        st["alk_alpha"] = results["alkalinity_alpha"] * 1e-6
        st["alk_beta"] = results["alkalinity_beta"] * 1e-6
        st["alkalinity_estimate"] = (
            st.alkalinity_from_pH
            + st.titrant_mass
            * self.titrant_molinity
            / (st.titrant_mass + self.analyte_mass)
        ) / st.dilution_factor
        if self.dic > 0:
            dic_loss = pyco2.sys(
                par1=st.alkalinity_from_acid.to_numpy() * 1e6,
                par2=st.pH.to_numpy(),
                par1_type=1,
                par2_type=3,
                opt_pH_scale=3,
                salinity=self.salinity,
                temperature=st.temperature.to_numpy(),
                **totals,
                **k_constants,
                uncertainty_from={"par1": 2, "par2": 0.01},
                uncertainty_into=["dic", "fCO2"],
            )
            st["dic_loss"] = dic_loss["dic"]
            st["dic_loss_u"] = dic_loss["u_dic"]
            st["dic_loss_lo"] = st.dic_loss - st.dic_loss_u
            st.dic_loss_lo.where(st.dic_loss_lo > 0, other=0, inplace=True)
            st["dic_loss_hi"] = st.dic_loss + st.dic_loss_u
            st["fCO2_loss"] = dic_loss["fCO2"]
            st["fCO2_loss_u"] = dic_loss["u_fCO2"]
            st["fCO2_loss_lo"] = st.fCO2_loss - st.fCO2_loss_u
            st.fCO2_loss_lo.where(st.fCO2_loss_lo > 0, other=0, inplace=True)
            st["fCO2_loss_hi"] = st.fCO2_loss + st.fCO2_loss_u
        else:
            st["dic_loss"] = 0.0
            st["dic_loss_u"] = 0.0
            st["dic_loss_lo"] = 0.0
            st["dic_loss_hi"] = 0.0
            st["fCO2_loss"] = 0.0
            st["fCO2_loss_u"] = 0.0
            st["fCO2_loss_lo"] = 0.0
            st["fCO2_loss_hi"] = 0.0

    def calkulate(
        self,
        alkalinity_certified,
        analyte_total_sulfate=None,
        emf0_guess=None,
        least_squares_kwargs=default.least_squares_kwargs,
        pH_range=None,
        titrant_molinity_guess=None,
        titrant=None,
    ):
        self.calibrate(
            alkalinity_certified,
            analyte_total_sulfate=analyte_total_sulfate,
            emf0_guess=emf0_guess,
            least_squares_kwargs=least_squares_kwargs,
            pH_range=pH_range,
            titrant_molinity_guess=titrant_molinity_guess,
            titrant=titrant,
        )
        self.solve(
            analyte_total_sulfate=analyte_total_sulfate,
            emf0_guess=emf0_guess,
            least_squares_kwargs=least_squares_kwargs,
            pH_range=pH_range,
            titrant=titrant,
        )

    def _get_dic_loss_hires(self, fCO2_air=None, split_pH=None):
        """Fit and forecast high-resolution DIC loss model."""
        self.check_set_default("fCO2_air", fCO2_air)
        self.check_set_default("split_pH", split_pH)
        k_dic_loss, loss_hires = core.loss.get_dic_loss_hires(
            self.titration.titrant_mass.to_numpy(),
            self.titration.pH.to_numpy(),
            self.titration.dic_loss.to_numpy(),
            self.titration.fCO2_loss.to_numpy(),
            self.titration.k_CO2.to_numpy(),
            self.titration.k_carbonic_1.to_numpy(),
            self.titration.k_carbonic_2.to_numpy(),
            self.analyte_mass,
            self.titration.dic.iloc[0],
            fCO2_air=self.fCO2_air,
            split_pH=self.split_pH,
        )
        return k_dic_loss, pd.DataFrame(loss_hires)

    def get_dic_loss(self, fCO2_air=None, split_pH=None):
        """Get final DIC loss values at the titration points for the titration df."""
        self.check_set_default("fCO2_air", fCO2_air)
        self.check_set_default("split_pH", split_pH)
        (
            self.k_dic_loss,
            self.titration["dic_loss_modelled"],
            self.titration["fCO2_loss_modelled"],
            self.titration["dic_loss_fitted"],
        ) = core.loss.get_dic_loss(
            self.titration.titrant_mass.to_numpy(),
            self.titration.pH.to_numpy(),
            self.titration.dic_loss.to_numpy(),
            self.titration.fCO2_loss.to_numpy(),
            self.titration.k_CO2.to_numpy(),
            self.titration.k_carbonic_1.to_numpy(),
            self.titration.k_carbonic_2.to_numpy(),
            self.analyte_mass,
            self.titration.dic.iloc[0],
            fCO2_air=self.fCO2_air,
            split_pH=self.split_pH,
        )

    def update_dic_loss(self, fCO2_air=None, split_pH=None):
        """Model and fit DIC-loss function and then re-solve accounting for DIC loss."""
        self.check_set_default("fCO2_air", fCO2_air)
        self.check_set_default("split_pH", split_pH)
        self.get_dic_loss()
        self.titration["dic"] = self.titration.dic_loss_modelled
        self.solve()

    # Assign plotting functions
    plot_emf = plot.titration.emf
    plot_pH = plot.titration.pH
    plot_gran_emf0 = plot.titration.gran_emf0
    plot_gran_alkalinity = plot.titration.gran_alkalinity
    plot_alkalinity = plot.titration.alkalinity
    plot_components = plot.titration.components
    plot_dic_loss = plot.titration.dic_loss
    plot_fCO2_loss = plot.titration.fCO2_loss

    def __str__(self):
        if self.file_name is not None:
            return f"Titration {self.file_name}"
        else:
            return "Titration (simulated)"

    def __repr__(self):
        rstr = "calkulate.Titration("
        if hasattr(self, "analyte_mass"):
            if self.analyte_mass is not None:
                rstr += f"\n    analyte_mass={self.analyte_mass},"
        if hasattr(self, "analyte_volume"):
            if self.analyte_volume is not None:
                rstr += f"\n    analyte_volume={self.analyte_volume},"
        if self.file_name is not None:
            rstr += f"\n    file_name='{self.file_name}',"
        if self.file_path is not None:
            rstr += f"\n    file_path='{self.file_path}',"
        if self.salinity != 35:
            rstr += f"\n    salinity={self.salinity},"
        rstr += "\n    **prepare_kwargs,\n)"
        return rstr
