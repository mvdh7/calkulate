<script type="text/x-mathjax-config">
MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}});
MathJax.Ajax.config.path["mhchem"] =
  "https://cdnjs.cloudflare.com/ajax/libs/mathjax-mhchem/3.3.2";
MathJax.Hub.Config({TeX: {extensions: ["[mhchem]/mhchem.js"]}});
</script><script src='https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-MML-AM_CHTML' async></script>

# Solving for alkalinity

Several different approaches for determining total alkalinity from titration data are implemented in Calkulate. They are implemented by functions in Calkulate's `solve` module, and they all share a common syntax:

```python
OptimizeResult = calk.solve.method(massAcid, emf, tempK, massSample, concAcid,
    concTotals, eqConstants)
```

The formats of the input variables are described in the documentation section on [variables and conventions](../conventions).

The output `OptimizeResult` is usually (except for the `halfGran` method, see later) the output of SciPy's `optimize.least_squares` function. As described in [the SciPy documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html), the actual fitted alkalinity value (and any auxiliary fitted parameters, which vary by method) is contained within the field `OptimizeResult['x']`. The alkalinity is always the first item in this list, i.e. `OptimizeResult['x'][0]`.

---

## Least-squares fitting

Modern alkalinity-solving methods typically involve setting up the alkalinity equation in terms of stoichiometric equilibrium constants and the hydrogen ion concentration, the latter of which can be calculated from the EMF measurements, given the electrode's EMF°. One then deploys a solver (e.g. least squares) to find the EMF° and alkalinity values that best fit the titration data points. The methods that follow this approach differ from each other in three main ways:

  1. Which chemical species are included in the alkalinity equation?
  2. Which titration points (i.e. what pH range) are used for the fitting?
  3. How is EMF° defined mathematically within the fitting process?

The different least-squares methods can thus be briefly defined in terms of the answers to the above questions, as follows.

### `complete`: complete calculation

This is the recommended method. It takes into account every major component of typical open ocean seawater, and uses a low pH range that minimises the influence of the carbonic acid equilibria upon the final result.

**Syntax:**

```python
alk, emf0 = calk.solve.complete(massAcid, emf, tempK, massSample, concAcid,
    concTotals, eqConstants)['x']
```

**Specifics:**

*1. Which chemical species are included in the alkalinity equation?*

> All of them! That is: $\ce{CO2(aq)}$, $\ce{HCO3-}$, $\ce{CO3^2-}$, $\ce{B(OH)4-}$, $\ce{HSO4-}$, $\ce{HF}$, $\ce{H3PO4}$, $\ce{HPO4^2-}$, $\ce{PO4^3-}$, $\ce{SiO(OH)3-}$, $\ce{H+}$ and $\ce{OH-}$.

*2. Which titration points (i.e. what pH range) are used for the fitting?*

> 3 < pH < 4

*3. How is EMF° defined mathematically within the fitting process?*

> Directly as EMF°.

### `DAA03`: Dickson CRM method

This is the approach described by [DAA03](../references/#DAA03). A separate implementation of this method is used to calibrate the widely used Certified Reference Material (CRM) seawater produced by Prof Andrew Dickson (Scripps Institution of Oceanography). It uses a low pH range and only includes chemical species that have appreciable concentrations in that range. Carbonic acid is ignored, as $\ce{CO2}$ is assumed to have been actively degassed before measurement.

**Syntax:**

```python
alk, f = calk.solve.DAA03(massAcid, emf, tempK, massSample, concAcid,
    concTotals, eqConstants)['x']
```

**Specifics:**

*1. Which chemical species are included in the alkalinity equation?*

> $\ce{HSO4-}$, $\ce{HF}$ and $\ce{H+}$.

*2. Which titration points (i.e. what pH range) are used for the fitting?*

> 3 < pH < 3.5

*3. How is EMF° defined mathematically within the fitting process?*

> As $f = [\ce{H+}]/[\ce{H'}]$, where $[\ce{H'}] = \exp[(\text{EMF°} - \text{EMF}) F / RT]$.

### `Dickson1981`: closed-cell titrations

This is the approach described by [D81](../references/#D81). It uses a high pH range and solves for dissolved inorganic carbon as well as total alkalinity, thus assuming that no $\ce{CO2}$ is degassed during the titration.

**Syntax:**

```python
alk, totalCarbonate, f = calk.solve.Dickson1981(massAcid, emf, tempK, massSample,
    concAcid, concTotals, eqConstants)['x']
```

Note that if a `totalCarbonate` value is provided within the `concTotals` dict it should not cause any problems but it will be completely ignored by this solver.

**Specifics:**

*1. Which chemical species are included in the alkalinity equation?*

> $\ce{CO2(aq)}$, $\ce{HCO3-}$, $\ce{CO3^2-}$, $\ce{B(OH)4-}$, $\ce{HSO4-}$, $\ce{HF}$, $\ce{H3PO4}$, $\ce{HPO4^2-}$, $\ce{PO4^3-}$, $\ce{H+}$ and $\ce{OH-}$.

*2. Which titration points (i.e. what pH range) are used for the fitting?*

> 5 < pH

*3. How is EMF° defined mathematically within the fitting process?*

> As $f = [\ce{H+}]/[\ce{H'}]$, where $[\ce{H'}] = \exp[(\text{EMF°} - \text{EMF}) F / RT]$.

---

## Gran plots

Prior to the advent of easy-to-compute least-squares solvers, alkalinity was determined from titration data using graphical approaches such as the modified Gran plot [[G52](../references/#G52); [HJ73](../references/#HJ73); [BBSW81](../references/#BBSW81)], which typically solve for not only alkalinity and EMF° but also dissolved inorganic carbon and the first dissociation constant of carbonic acid. Only the 'half-Gran' method of [H15](../references/#H15) is implemented in Calkulate, which instead takes dissolved inorganic carbon as an input and solves only for alkalinity and EMF°.

### `halfGran`: half-Gran method

**Syntax:**

```python
alk, emf0 = calk.solve.halfGran(massAcid, emf, tempK, massSample, concAcid,
    concTotals, eqConstants)['x']
```

Note that this returns a dict with one field (`'x'`) containing a list of the alkalinity and EMF° values, rather than a true `OptimizeResult` as in the case of the least-squares fitting functions.

**Specifics:**

*1. Which chemical species are included in the alkalinity equation?*

> $\ce{CO2(aq)}$, $\ce{HCO3-}$, $\ce{CO3^2-}$, $\ce{B(OH)4-}$, $\ce{HSO4-}$, $\ce{HF}$, $\ce{H3PO4}$, $\ce{HPO4^2-}$, $\ce{PO4^3-}$, $\ce{H+}$ and $\ce{OH-}$.

*2. Which titration points (i.e. what pH range) are used for the fitting?*

> 3 < pH < 4

*3. How is EMF° defined mathematically within the fitting process?*

> Directly as EMF°.
