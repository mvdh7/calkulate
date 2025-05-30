import os

import pandas as pd

from .dataset import _backcompat, _get_kwargs_for, calibrate_row, solve_row
from .read.titrations import keys_read_dat, read_dat
from .titration import (
    convert_amount_units,
    get_totals_k_constants,
    keys_cau,
    keys_totals_ks,
)


class Dataset(pd.DataFrame):
    """
    `calk.Dataset`
    ==============
    A `calk.Dataset` is pandas `DataFrame` with several `calk.dataset`
    functions available as methods.

    Alkalinity solving methods
    --------------------------
    calibrate
        Find the best-fit `titrant_molinity` for each sample that has a value
        for `alkalinity_certified`.
    solve
        Solve every sample for `alkalinity` when `titrant_molinity` is known.
    calkulate
        Run the `calibrate` and `solve` steps sequentially.

    Data visualisation methods
    --------------------------
    plot_titrant_molinity
        Plot the `titrant_molinity` values determined with `calibrate` through time.
    plot_alkalinity_offset
        Plot the difference between solved `alkalinity` and `alkalinity_certified`
        through time.

    Conversion methods
    ------------------
    to_Titration
        Return a `calk.Titration` object for one row in the `Dataset`.
    to_pandas
        Return a copy of the `Dataset` as a standard pandas `DataFrame`.
    """

    from .dataset import calibrate, calkulate, solve

    # get_batches = get_batches
    # get_total_salts = get_total_salts
    # to_Titration = to_Titration

    # from .plot import (
    #     alkalinity_offset as plot_alkalinity_offset,
    #     titrant_molinity as plot_titrant_molinity,
    # )

    def to_pandas(self):
        """Return a copy of the `Dataset` as a standard pandas `DataFrame`."""
        return pd.DataFrame(self)


class Titration:
    """
    `calk.Titration`
    ================
    A structure containing all the information about one titration with methods for
    calibrating and solving and then visualising the results.

    Creating a `Titration`
    ----------------------
    There are two ways to create a `Titration`:

    (1) From a row in `calk.Dataset` using the `to_Titration` method (or
    equivalently from a `pandas.DataFrame` using the `calk.to_Titration`
    function):

      >>> tt = ds.to_Titration(index)
      >>> tt = calk.to_Titration(ds, index)

    (2) By initalising a `calk.Titration` directly with a `file_name` for a
    titration data file.

      >>> tt = calk.Titration(file_name="path/to/file.dat", **kwargs)

    Alkalinity solving methods
    --------------------------
    `calibrate`
        Find the best-fit `titrant_molinity` for a given `alkalinity_certified`.
    `set_titrant_molinity`
        Set the `titrant_molinity` if it is known independently.
    `solve`
        Solve for `alkalinity` when `titrant_molinity` is known.
    `calkulate`
        Run the `calibrate` and `solve` steps sequentially.

    Data visualisation methods
    --------------------------
    `plot_emf`
        EMF throughout the titration.
    `plot_pH`
        pH throughout the titration.
    `plot_gran_emf0`
        The Gran-plot initial estimate of EMF0.
    `plot_gran_alkalinity`
        The Gran-plot initial estimate of alkalinity.
    `plot_alkalinity`
        The final alkalinity estimates across the titration.
    `plot_components`
        Show how all equilibrating components of the solution change throughout the
        titration.

    Attributes
    ----------
    `alkalinity` : `float`
        The solved alkalinity in micromol/kg-seawater.
    `alkalinity_certified` : `float`
        The certified alkalinity value against which the sample was calibrated,
        if applicable, in micromol/kg-seawater.
    `alkalinity_gran` : `float`
        The Gran-plot initial estimate of alkalinity in micromol/kg-seawater.
    `analyte_mass` : `float`
        The mass of the analyte in kg.
    `analyte_volume` : `float`
        The volume of the analyte in ml, if provided.
    `calibrated` : `bool`
        Whether or not the `Titration` has been calibrated.
    `dic` : `float`
        The dissolved inorganic carbon in micromol/kg-seawater.
    `emf0` : `float`
        The solved EMF0 in mV.
    `emf0_gran` : `float`
        The Gran-plot initial estimate of EMF0 in mV.
    `file_name` : `str`
        The name of the titration data file.
    `file_path` : `str`
        The path to the titration data file.
    `file_prepare_kwargs` : `dict`
        The kwargs passed to prepare the titration data file for solving.
    `opt_result` : `scipy.optimize._optimize.OptimizeResult`
        The raw output from the `scipy.optimize.least_squares` solver.
    `pH_initial` : `float`
        The pH at the initial titration point on the free pH scale.
    `pH_initial_temperature` : `float`
        The temperature at the initial titration point.
    `pH_range` : `tuple` or `list`
        The (minimum, maximum) pH to use for solving the titration.
    `salinity` : `float`
        The practical salinity of the analyte.
    `solved` : `bool`
        Whether or not the titration has been solved.
    `titrant` : `str`
        The type of titrant used in the titration.
    `titrant_molinity` : `float`
        The molinity of the titrant in mol/kg-titrant.
    `titration` : `pd.DataFrame`
        A table of everything that changes through the titration (e.g., amount of
        titrant, EMF, temperature, total salt concentrations and equilibrium constants).
    """

    def __init__(self, row, **kwargs):
        self.row = row.copy()
        self.kwargs = _backcompat(row, kwargs.copy())
        self.titration = {}
        if "file_good" not in row:
            self.row["file_good"] = True
        self.get_titration()
        if "titrant_molinity" in self.row and pd.notnull(
            self.row.titrant_molinity
        ):
            self.solve()
        elif "alkalinity_certified" in self.row and pd.notnull(
            self.row.alkalinity_certified
        ):
            self.calibrate().solve()

    def __getattr__(self, attr):
        # Enable using self.attr as a shortcut for self.row.attr
        try:
            return object.__getattribute__(self, attr)
        except AttributeError:
            try:
                return self.row[attr]
            except KeyError:
                raise AttributeError(attr)

    def __getitem__(self, key):
        # Enable using self["attr"] as a shortcut for self.row["attr"]
        try:
            return self.row[key]
        except KeyError:
            raise AttributeError(key)

    def get_titration(self):
        # STEP 1: IMPORT TITRATION DATA
        # Append file_name to file_path if present
        if "file_path" in self.row:
            file_name = os.path.join(self.row.file_path, self.row.file_name)
        else:
            file_name = self.row.file_name
        # Get kwargs for read_dat
        kwargs_read_dat = _get_kwargs_for(keys_read_dat, self.kwargs, self.row)
        # Import the file
        try:
            dat_data = read_dat(file_name, **kwargs_read_dat)
        except FileNotFoundError:
            raise FileNotFoundError(f'File not found: "{file_name}"')
        for k, v in dat_data._asdict().items():
            self.titration[k] = v
        # STEP 2: CONVERT AMOUNT UNITS
        # Get kwargs for convert_amount_units
        kwargs_cau = _get_kwargs_for(keys_cau, self.kwargs, self.row)
        # If there is analyte_mass, use it, otherwise use analyte_volume
        if "analyte_volume" in kwargs_cau and "analyte_mass" in kwargs_cau:
            kwargs_cau["analyte_volume"] = None
        # Prepare titration file, totals and k_constants
        try:
            converted = convert_amount_units(
                dat_data, self.row.salinity, **kwargs_cau
            )
            self.row["analyte_mass"] = converted.analyte_mass
        except Exception as e:
            print(f'Error converting units for "{file_name}":')
            print(f"{e}")
        for k, v in converted._asdict().items():
            self.titration[k] = v
        # STEP 3: GET TOTAL SALTS AND EQUILIBRIUM CONSTANTS
        kwargs_totals_ks = _get_kwargs_for(
            keys_totals_ks, self.kwargs, self.row
        )
        try:
            totals, k_constants = get_totals_k_constants(
                converted, self.row.salinity, **kwargs_totals_ks
            )
        except Exception as e:
            print(f'Error getting salts and constants for "{file_name}":')
            print(f"{e}")
        for k, v in totals.items():
            self.titration[k] = v
        for k, v in k_constants.items():
            self.titration[k] = v
        self.titration = pd.DataFrame(self.titration)
        return self

    def calibrate(self):
        calibrated = calibrate_row(self.row, **self.kwargs)
        for k, v in calibrated.items():
            self.row[k] = v
        self.row["titrant_molinity"] = self.row["titrant_molinity_here"]
        return self

    def solve(self):
        solved = solve_row(self.row, **self.kwargs)
        for k, v in solved.items():
            self.row[k] = v
        return self

    # @np.errstate(invalid="ignore")
    # def gran_guesses(self, emf0_guess=None):
    #     # Get simple Gran-plot estimator
    #     st = self.titration
    #     st["mixture_mass"] = st.titrant_mass + self.analyte_mass
    #     st["gran_estimates"] = core.gran_estimator(
    #         st.mixture_mass, st.emf, st.temperature
    #     )
    #     # Select which data points to use for first guesses
    #     st["G_gran"] = (
    #         st.gran_estimates >= 0.1 * np.max(st.gran_estimates)
    #     ) & (st.gran_estimates <= 0.9 * np.max(st.gran_estimates))
    #     # Make first guesses
    #     (
    #         alkalinity_gran,
    #         self.gran_slope,
    #         self.gran_intercept,
    #     ) = core.gran_guess_alkalinity(
    #         st.titrant_mass[st.G_gran],
    #         st.gran_estimates[st.G_gran],
    #         self.analyte_mass,
    #         self.titrant_molinity,
    #     )
    #     if self.titrant == "H2SO4":
    #         alkalinity_gran *= 2
    #     if emf0_guess is None:
    #         st["emf0_gran"] = core.gran_guesses_emf0(
    #             st.titrant_mass,
    #             st.emf,
    #             st.temperature,
    #             self.analyte_mass,
    #             self.titrant_molinity,
    #             alkalinity_guess=alkalinity_gran,
    #             titrant=self.titrant,
    #             HF=0,
    #             HSO4=0,
    #         )
    #         self.emf0_gran = np.mean(st.emf0_gran[st.G_gran])
    #     else:
    #         self.emf0_gran = emf0_guess
    #     st["pH_gran"] = convert.emf_to_pH(
    #         st.emf, self.emf0_gran, st.temperature
    #     )
    #     self.alkalinity_gran = alkalinity_gran * 1e6
    #     self.titration["G_final"] = (
    #         self.titration.pH_gran >= self.pH_range[0]
    #     ) & (self.titration.pH_gran <= self.pH_range[1])

    # def solve(
    #     self,
    #     analyte_total_sulfate=None,
    #     emf0_guess=None,
    #     least_squares_kwargs=default.least_squares_kwargs,
    #     pH_range=None,
    #     titrant_molinity=None,
    #     titrant=None,
    # ):
    #     if analyte_total_sulfate is not None:
    #         self.analyte_total_sulfate = analyte_total_sulfate
    #     if titrant_molinity is not None:
    #         self.set_titrant_molinity(titrant_molinity)
    #     self.check_set_default("pH_range", pH_range)
    #     self.check_set_default("titrant", titrant)
    #     self.gran_guesses()
    #     # Solve for total alkalinity and EMF0
    #     st = self.titration
    #     if self.titrant == "H2SO4":
    #         # Update sulfate dilution by titrant and its consequences
    #         assert self.analyte_total_sulfate is not None
    #         st["total_sulfate"] = (
    #             1e-6 * self.analyte_total_sulfate * self.analyte_mass
    #             + st.titrant_molinity * st.titrant_mass
    #         ) / (self.analyte_mass + st.titrant_mass)
    #         # Solve for alkalinity and EMF0
    #         self.opt_result = core.solve_emf_complete_H2SO4(
    #             self.titrant_molinity,
    #             st.titrant_mass.to_numpy(),
    #             st.emf.to_numpy(),
    #             st.temperature.to_numpy(),
    #             self.analyte_mass,
    #             self.get_totals(),
    #             self.get_k_constants(),
    #             least_squares_kwargs=least_squares_kwargs,
    #             pH_range=self.pH_range,
    #             emf0_guess=emf0_guess,
    #         )
    #     else:
    #         self.opt_result = core.solve_emf_complete(
    #             self.titrant_molinity,
    #             st.titrant_mass.to_numpy(),
    #             st.emf.to_numpy(),
    #             st.temperature.to_numpy(),
    #             self.analyte_mass,
    #             self.get_totals(),
    #             self.get_k_constants(),
    #             least_squares_kwargs=least_squares_kwargs,
    #             pH_range=self.pH_range,
    #             emf0_guess=emf0_guess,
    #         )
    #     self.alkalinity, self.emf0 = self.opt_result["x"]
    #     self.alkalinity *= 1e6
    #     # Calculate initial pH
    #     st["pH"] = convert.emf_to_pH(st.emf, self.emf0, st.temperature)
    #     self.pH_initial = st["pH"].iloc[0]
    #     self.pH_initial_temperature = st["temperature"].iloc[0]
    #     self.get_alkalinity_from_acid()
    #     self.do_CO2SYS()
    #     self.solved = True

    # def get_alkalinity_from_acid(self):
    #     self.titration["alkalinity_from_acid"] = (
    #         1e-6 * self.alkalinity * self.analyte_mass
    #         - self.titration.titrant_mass * self.titrant_molinity
    #     ) / (self.analyte_mass + self.titration.titrant_mass)

    # def do_CO2SYS(self):
    #     st = self.titration
    #     totals = {
    #         k: st[k].to_numpy() * 1e6 if k in st else 0
    #         for k in [
    #             "total_sulfate",
    #             "total_fluoride",
    #             "total_borate",
    #             "total_phosphate",
    #             "total_silicate",
    #             "total_ammonia",
    #             "total_sulfide",
    #             "total_alpha",
    #             "total_beta",
    #         ]
    #     }
    #     k_constants = {
    #         k: st[k].to_numpy() for k in st.columns if k.startswith("k_")
    #     }
    #     results = pyco2.sys(
    #         par1=st.dic.to_numpy() * 1e6,
    #         par2=st.pH.to_numpy(),
    #         par1_type=2,
    #         par2_type=3,
    #         opt_pH_scale=3,
    #         salinity=self.salinity,
    #         temperature=st.temperature.to_numpy(),
    #         **totals,
    #         **k_constants,
    #     )
    #     st["alkalinity_from_pH"] = results["alkalinity"] * 1e-6
    #     st["k_CO2"] = results["k_CO2"]
    #     for co2sysvar in [
    #         "HCO3",
    #         "CO3",
    #         "BOH4",
    #         "PO4",
    #         "HPO4",
    #         "H3PO4",
    #         "HSO4",
    #         "HF",
    #         "H3SiO4",
    #         "NH3",
    #         "HS",
    #         "OH",
    #     ]:
    #         st[co2sysvar] = results[co2sysvar] * 1e-6
    #     st["H"] = 10**-st.pH
    #     st["alk_alpha"] = results["alkalinity_alpha"] * 1e-6
    #     st["alk_beta"] = results["alkalinity_beta"] * 1e-6
    #     st["alkalinity_estimate"] = (
    #         st.alkalinity_from_pH
    #         + st.titrant_mass
    #         * self.titrant_molinity
    #         / (st.titrant_mass + self.analyte_mass)
    #     ) / st.dilution_factor
    #     if self.dic > 0:
    #         dic_loss = pyco2.sys(
    #             par1=st.alkalinity_from_acid.to_numpy() * 1e6,
    #             par2=st.pH.to_numpy(),
    #             par1_type=1,
    #             par2_type=3,
    #             opt_pH_scale=3,
    #             salinity=self.salinity,
    #             temperature=st.temperature.to_numpy(),
    #             **totals,
    #             **k_constants,
    #             uncertainty_from={"par1": 2, "par2": 0.01},
    #             uncertainty_into=["dic", "fCO2"],
    #         )
    #         st["dic_loss"] = dic_loss["dic"]
    #         st["dic_loss_u"] = dic_loss["u_dic"]
    #         st["dic_loss_lo"] = st.dic_loss - st.dic_loss_u
    #         st.dic_loss_lo.where(st.dic_loss_lo > 0, other=0, inplace=True)
    #         st["dic_loss_hi"] = st.dic_loss + st.dic_loss_u
    #         st["fCO2_loss"] = dic_loss["fCO2"]
    #         st["fCO2_loss_u"] = dic_loss["u_fCO2"]
    #         st["fCO2_loss_lo"] = st.fCO2_loss - st.fCO2_loss_u
    #         st.fCO2_loss_lo.where(st.fCO2_loss_lo > 0, other=0, inplace=True)
    #         st["fCO2_loss_hi"] = st.fCO2_loss + st.fCO2_loss_u
    #     else:
    #         st["dic_loss"] = 0.0
    #         st["dic_loss_u"] = 0.0
    #         st["dic_loss_lo"] = 0.0
    #         st["dic_loss_hi"] = 0.0
    #         st["fCO2_loss"] = 0.0
    #         st["fCO2_loss_u"] = 0.0
    #         st["fCO2_loss_lo"] = 0.0
    #         st["fCO2_loss_hi"] = 0.0

    # def _get_dic_loss_hires(self, fCO2_air=None, split_pH=None):
    #     """Fit and forecast high-resolution DIC loss model."""
    #     self.check_set_default("fCO2_air", fCO2_air)
    #     self.check_set_default("split_pH", split_pH)
    #     k_dic_loss, loss_hires = core.loss.get_dic_loss_hires(
    #         self.titration.titrant_mass.to_numpy(),
    #         self.titration.pH.to_numpy(),
    #         self.titration.dic_loss.to_numpy(),
    #         self.titration.fCO2_loss.to_numpy(),
    #         self.titration.k_CO2.to_numpy(),
    #         self.titration.k_carbonic_1.to_numpy(),
    #         self.titration.k_carbonic_2.to_numpy(),
    #         self.analyte_mass,
    #         self.titration.dic.iloc[0],
    #         fCO2_air=self.fCO2_air,
    #         split_pH=self.split_pH,
    #     )
    #     return k_dic_loss, pd.DataFrame(loss_hires)

    # def get_dic_loss(self, fCO2_air=None, split_pH=None):
    #     """Get final DIC loss values at the titration points for the titration df."""
    #     self.check_set_default("fCO2_air", fCO2_air)
    #     self.check_set_default("split_pH", split_pH)
    #     (
    #         self.k_dic_loss,
    #         self.titration["dic_loss_modelled"],
    #         self.titration["fCO2_loss_modelled"],
    #         self.titration["dic_loss_fitted"],
    #     ) = core.loss.get_dic_loss(
    #         self.titration.titrant_mass.to_numpy(),
    #         self.titration.pH.to_numpy(),
    #         self.titration.dic_loss.to_numpy(),
    #         self.titration.fCO2_loss.to_numpy(),
    #         self.titration.k_CO2.to_numpy(),
    #         self.titration.k_carbonic_1.to_numpy(),
    #         self.titration.k_carbonic_2.to_numpy(),
    #         self.analyte_mass,
    #         self.titration.dic.iloc[0],
    #         fCO2_air=self.fCO2_air,
    #         split_pH=self.split_pH,
    #     )

    # def update_dic_loss(self, fCO2_air=None, split_pH=None):
    #     """Model and fit DIC-loss function and then re-solve accounting for DIC loss."""
    #     self.check_set_default("fCO2_air", fCO2_air)
    #     self.check_set_default("split_pH", split_pH)
    #     self.get_dic_loss()
    #     self.titration["dic"] = self.titration.dic_loss_modelled
    #     self.solve()

    # # Assign plotting functions
    # plot_emf = plot.titration.emf
    # plot_pH = plot.titration.pH
    # plot_gran_emf0 = plot.titration.gran_emf0
    # plot_gran_alkalinity = plot.titration.gran_alkalinity
    # plot_alkalinity = plot.titration.alkalinity
    # plot_components = plot.titration.components
    # plot_dic_loss = plot.titration.dic_loss
    # plot_fCO2_loss = plot.titration.fCO2_loss

    # def __str__(self):
    #     if self.file_name is not None:
    #         return f"Titration {self.file_name}"
    #     else:
    #         return "Titration (simulated)"

    # def __repr__(self):
    #     rstr = "calkulate.Titration("
    #     if hasattr(self, "analyte_mass"):
    #         if self.analyte_mass is not None:
    #             rstr += f"\n    analyte_mass={self.analyte_mass},"
    #     if hasattr(self, "analyte_volume"):
    #         if self.analyte_volume is not None:
    #             rstr += f"\n    analyte_volume={self.analyte_volume},"
    #     if self.file_name is not None:
    #         rstr += f"\n    file_name='{self.file_name}',"
    #     if self.file_path is not None:
    #         rstr += f"\n    file_path='{self.file_path}',"
    #     if self.salinity != 35:
    #         rstr += f"\n    salinity={self.salinity},"
    #     rstr += "\n    **prepare_kwargs,\n)"
    #     return rstr
