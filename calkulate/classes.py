import os

import pandas as pd
import PyCO2SYS as pyco2

from .convert import get_dilution_factor
from .core import SolveEmfResult
from .dataset import _backcompat
from .files import solve


def to_Titration(ds, index, **kwargs):
    """Create a `calk.Titration` for one titration in a dataset.

    Parameters
    ----------
    ds : pandas.DataFrame or calk.Dataset
        The dataset to make the Titration from.
    index
        The row index in the `ds` to use.
    kwargs
        Any additional kwargs, which will overwrite equivalent keys in the
        `ds` row.

    Returns
    -------
    calk.Titration
        A `calk.Titration` for the specified row of the `Dataset`.
    """
    return Titration(ds.loc[index], **kwargs)


class Dataset(pd.DataFrame):
    """
    calkulate.Dataset
    =================
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
        Plot the `titrant_molinity` values determined with `calibrate` through
        time.
    plot_alkalinity_offset
        Plot the difference between solved `alkalinity` and
        `alkalinity_certified` through time.

    Conversion methods
    ------------------
    to_Titration
        Return a `calk.Titration` for one titration in the `Dataset`.
    to_pandas
        Return a copy of the `Dataset` as a pandas `DataFrame`.
    """

    from .dataset import calibrate, calkulate, solve

    def to_Titration(self, index, **kwargs):
        """Create a `calk.Titration` for one titration in the dataset.

        Parameters
        ----------
        index
            The row index in the dataset to use.
        kwargs
            Any additional kwargs, which will overwrite equivalent keys in the
            dataset row.

        Returns
        -------
        calk.Titration
            A `calk.Titration` for the specified row of the dataset.
        """
        return to_Titration(self, index, **kwargs)

    from .plot import (
        alkalinity_offset as plot_alkalinity_offset,
        titrant_molinity as plot_titrant_molinity,
    )

    def to_pandas(self):
        """Return a copy of the `Dataset` as a standard pandas `DataFrame`."""
        return pd.DataFrame(self)


class Titration:
    """
    calkulate.Titration
    ===================
    A structure containing all the information about one titration with methods
    for calibrating and solving and then visualising the results.

    Creating a `Titration`
    ----------------------
    There are two ways to create a `Titration`:

    (1) From a row in `calk.Dataset` using the `to_Titration` method (or
    equivalently from a `pandas.DataFrame` using the `calk.to_Titration`
    function):

      >>> index = 0
      >>> tt = ds.to_Titration(index)
      >>> tt = calk.to_Titration(ds, index)

    (2) By initalising a `calk.Titration` directly with a dataset row:

      >>> tt = calk.Titration(ds.loc[index], **kwargs)

    Data visualisation methods
    --------------------------
    plot_emf
        EMF throughout the titration.
    plot_pH
        pH throughout the titration.
    plot_gran_emf0
        The Gran-plot initial estimate of EMF0.
    plot_gran_alkalinity
        The Gran-plot initial estimate of alkalinity.
    plot_alkalinity
        The final alkalinity estimates across the titration.
    plot_components
        Show how all equilibrating components of the solution change throughout the
        titration.

    Attributes
    ----------
    alkalinity : float
        The solved alkalinity in µmol/kg-seawater.
    alkalinity_certified : float
        The certified alkalinity value against which the sample was calibrated,
        if applicable, in µmol/kg-seawater.
    analyte_mass : float
        The mass of the analyte in kg.
    analyte_volume : float
        The volume of the analyte in ml, if provided.
    dic : float
        The dissolved inorganic carbon in µmol/kg-seawater.
    emf0 : float
        The solved EMF0 in mV.
    file_name : str
        The name of the titration data file.
    file_path : str
        The path to the titration data file.
    gran_alkalinity : float
        The Gran-plot initial estimate of alkalinity in µmol/kg-seawater.
    gran_emf0 : float
        The Gran-plot or user-provided initial estimate of EMF0 in mV.
    opt_result : scipy.optimize._optimize.OptimizeResult
        The raw output from the `scipy.optimize.least_squares` solver.
    salinity : float
        The practical salinity of the analyte.
    titrant_molinity : float
        The molinity of the titrant in mol/kg-titrant.
    titration : pd.DataFrame
        A table of everything that changes through the titration (e.g., amount
        of titrant, EMF, temperature, total salt concentrations and equilibrium
        constants).  `t` can be used as a shortcut for this.
    """

    def __init__(self, row, **kwargs):
        self.row = row.copy()
        self.kwargs = _backcompat(kwargs, row)
        for k, v in kwargs.items():
            self.row[k] = v
        self.sr = solve(**row)
        self.t = {
            "titrant_mass": self.sr.titrant_mass,
            "temperature": self.sr.temperature,
            "pH": self.sr.pH,
            "used": self.sr.used,
        }
        self.row["analyte_mass"] = self.sr.analyte_mass
        self.row["alkalinity"] = self.sr.alkalinity
        self.row["alkalinity_npts"] = self.sr.used.sum()
        self.row["alkalinity_std"] = self.sr.alkalinity_all[self.sr.used].std()
        self.row["pH_init"] = self.sr.pH[0]
        self.row["temperature_init"] = self.sr.temperature[0]
        if isinstance(self.sr, SolveEmfResult):
            self.row["gran_alkalinity"] = self.sr.ggr.alkalinity * 1e6
            self.row["gran_emf0"] = self.sr.ggr.emf0
            self.row["emf0"] = self.sr.emf0
            self.t.update(
                {
                    "emf": self.sr.emf,
                    "gran_func": self.sr.ggr.gfunc,
                    "gran_emf0": self.sr.ggr.emf0s,
                    "gran_pH": self.sr.ggr.pH,
                    "gran_used": self.sr.ggr.used,
                }
            )
        self.t.update({**self.sr.totals, **self.sr.k_constants})
        self.t = self.titration = pd.DataFrame(self.t)
        self.get_alkalinity_from_acid()
        self.do_CO2SYS()

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

    def __str__(self):
        if "file_path" in self.row:
            f = os.path.join(self.row.file_path, self.row.file_name)
        else:
            f = self.ro.file_name
        return f'calkulate.Titration for file "{f}"'

    def __repr__(self):
        if "file_path" in self.row:
            f = os.path.join(self.row.file_path, self.row.file_name)
        else:
            f = self.ro.file_name
        return f'calkulate.Titration for file "{f}"'

    def get_alkalinity_from_acid(self):
        self.titration["alkalinity_from_acid"] = (
            1e-6 * self.alkalinity * self.analyte_mass
            - self.titration.titrant_mass * self.titrant_molinity
        ) / (self.analyte_mass + self.titration.titrant_mass)

    def do_CO2SYS(self):
        st = self.titration
        st["dilution_factor"] = get_dilution_factor(
            st.titrant_mass, self.row.analyte_mass
        )
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
        k_constants = {
            k: st[k].to_numpy() for k in st.columns if k.startswith("k_")
        }
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
            "NH3",
            "HS",
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

    from .plot.titration import (
        alkalinity as plot_alkalinity,
        components as plot_components,
        emf as plot_emf,
        gran_alkalinity as plot_alkalinity_gran,  # TODO add trend line
        gran_emf0 as plot_emf0_gran,
        pH as plot_pH,
    )
