# Data import and export

`.io` contains functions to import titration data and export results in various formats.

<hr />

## `.vindta` - import VINDTA-style .dat files

Imports titration data in the .dat file format output by VINDTA instruments. This format is easy to replicate with any other sort of data that you may have.

**Syntax:**

```python
Vacid, EMF, tempK = calk.io.vindta(datfile)
```

**Input:**

  * `datfile` - string containing the .dat file's name (and path);

**Outputs:**

  * `Vacid` - volume of acid titrant added in ml;
  * `EMF` - titration cell potential in mV;
  * `tempK` - titre temperature in K.

<hr />

## `.Dickson1981` - load simulated titration

Imports data from titration simulated by Dickson (1981), reported in his Tables 1 (without phosphate) and 4 (with phosphate).

**Syntax:**

```python
Macid, pH, tempK, Msamp, Cacid, psal, XT, KX = calk.io.Dickson1981(withPhosphate=True)
```

**Input:**

  * `withPhosphate` - logical determining whether titration data should include phosphate (`True`, Table 4) or not (`False`, Table 1); defaults to `True`.

**Outputs:**

  * `Macid`, `tempK`, `Msamp`, `Cacid`, `psal`, `XT` and `KX` have [their usual meanings]();
  * `XT[0]` contains the true total alkalinity value;
  * `pH` - the pH values (Free scale) simulated by Dickson (1981) throughout the titration.
