# %%
import pandas as pd
import PyCO2SYS as pyco2

import calkulate as calk
from calkulate import files
from calkulate.convert import get_dilution_factor
from calkulate.core import SolveEmfResult
from calkulate.dataset import _backcompat


ds = calk.read_csv("tests/data/oberlin/metadf.csv")
row = ds.iloc[0]

# sr = calk.files.solve(**row)


class Titration:
    def __init__(self, row, **kwargs):
        self.row = row.copy()
        self.kwargs = _backcompat(kwargs, row)
        self.sr = files.solve(**row, **kwargs)
        self.t = {
            "titrant_mass": self.sr.titrant_mass,
            "temperature": self.sr.temperature,
            "pH": self.sr.pH,
            "used": self.sr.used,
        }
        self.row["alkalinity"] = self.sr.alkalinity
        self.row["emf0"] = self.sr.emf0
        self.row["analyte_mass"] = self.sr.analyte_mass
        if isinstance(self.sr, SolveEmfResult):
            self.row["gran_alkalinity"] = self.sr.ggr.alkalinity
            self.row["gran_emf0"] = self.sr.ggr.emf0
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

    from calkulate.plot.titration import (
        alkalinity as plot_alkalinity,
        components as plot_components,
        emf as plot_emf,
        gran_alkalinity as plot_alkalinity_gran,  # TODO add trend line
        gran_emf0 as plot_emf0_gran,
        pH as plot_pH,
    )


tt = calk.classes.Titration(row, dic=650)
tt.t
# tt.plot_emf()
# tt.plot_pH()
# tt.plot_emf0_gran()
tt.plot_alkalinity_gran()
# tt.plot_alkalinity()
# tt.plot_components()
