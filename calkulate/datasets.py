# Calkulate: seawater total alkalinity from titration data
# Copyright (C) 2019-2020  Matthew P. Humphreys  (GNU GPLv3)
"""Import and export titration data."""

import numpy as np
import pandas as pd
from . import types


class Dataset:
    """A full titration dataset."""

    def __init__(self, df, get_titrations=True, read_func=None, **read_kwargs):
        if isinstance(df, str):
            if read_func is not None:
                self.table = read_func(df, **read_kwargs)
            elif df.endswith(".csv"):
                self.table = pd.read_csv(df, **read_kwargs)
            elif df.endswith(".xlsx") or df.endswith(".xls"):
                self.table = pd.read_excel(df, **read_kwargs)
        else:
            self.table = pd.DataFrame(df)
        self.titrations = {}
        if get_titrations:
            self.import_titrations()
        if "analysis_batch" not in self.table:
            self.table["analysis_batch"] = 0
        self.batch_groups = self.table.groupby(by="analysis_batch")
        self.batches = self.batch_groups.analysis_batch.agg(analysis_count="count")

    def import_titrations(self):
        any_errors = False
        if "file_good" not in self.table.columns:
            self.table["file_good"] = True
        for i in self.table.index:
            if self.table.loc[i].file_good:
                ttr = self.table.loc[i]
                if "file_path" in ttr:
                    fname = ttr.file_path + ttr.file_name
                else:
                    fname = ttr.file_name
                try:
                    self.titrations.update({i: types.Titration(ttr)})
                except IOError:
                    print("Can't find file: '{}'.".format(fname))
                    self.titrations.update({i: None})
                    any_errors = True
                except:
                    print("Error importing file: '{}'.".format(fname))
                    self.titrations.update({i: None})
                    any_errors = True
            else:
                self.titrations.update({i: None})
        if any_errors:
            print("All other titration files imported successfully!")
        else:
            print("All titration files imported successfully!")

    def calibrate_titrants(self):
        """Calibrate all titrations that have a certified alkalinity value."""
        assert (
            "alkalinity_certified" in self.table
        ), "Missing 'alkalinity_certified' field."
        if "titrant_molinity_calibrated" not in self.table:
            self.table["titrant_molinity_calibrated"] = np.nan
        for i in self.table.index:
            if self.titrations[i] is not None:
                if ~np.isnan(self.table.loc[i].alkalinity_certified):
                    self.titrations[i].calibrate()
                    self.table.loc[i, "titrant_molinity_calibrated"] = self.titrations[
                        i
                    ].titrant.molinity_calibrated

    def calibrate_batches(self):
        """Assemble calibrated titrant molinities by batch and broadcast into table."""
        self.batches = self.batches.join(
            (
                self.batch_groups.titrant_molinity_calibrated.agg(
                    titrant_molinity=np.mean,
                    titrant_molinity_std=np.std,
                    titrant_molinity_count=lambda x: np.sum(~np.isnan(x)),
                ),
            )[0]
        )
        self.table["titrant_molinity"] = self.batches.loc[
            self.table.analysis_batch
        ].titrant_molinity.values

    def calibrate(self):
        """Perform all titrant calibration steps."""
        self.calibrate_titrants()
        self.calibrate_batches()

    def solve(self):
        """Solve all titrations for alkalinity."""
        assert "titrant_molinity" in self.table, "Missing 'titrant_molinity' field."
        if "alkalinity" not in self.table:
            self.table["alkalinity"] = np.nan
            self.table["emf0"] = np.nan
            self.table["pH"] = np.nan
            self.table["pH_temperature"] = np.nan
        for i in self.table.index:
            tti = self.titrations[i]
            if tti is not None:
                if ~np.isnan(self.table.loc[i].titrant_molinity):
                    try:
                        tti.titrant.molinity = self.table.loc[i].titrant_molinity
                        tti.solve()
                        self.table.loc[i, "alkalinity"] = tti.analyte.alkalinity
                        if tti.measurement_type == "EMF":
                            self.table.loc[i, "emf0"] = tti.analyte.emf0
                            self.table.loc[i, "pH"] = tti.analyte.pH
                            self.table.loc[
                                i, "pH_temperature"
                            ] = tti.analyte.pH_temperature
                    except:
                        print("Failed to solve: '{}'".format(tti.fname))

    def calibrate_and_solve(self):
        """Perform all titrant calibration steps and solve all titrations for alkalinity."""
        self.calibrate()
        self.solve()
