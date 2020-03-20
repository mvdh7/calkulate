<script type="text/x-mathjax-config">
MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}});
MathJax.Ajax.config.path["mhchem"] =
  "https://cdnjs.cloudflare.com/ajax/libs/mathjax-mhchem/3.3.2";
MathJax.Hub.Config({TeX: {extensions: ["[mhchem]/mhchem.js"]}});
</script><script src='https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-MML-AM_CHTML' async></script>

# Total solute concentrations

The [alkalinity solvers](../solvers) need to know the total concentrations of chemical species involved in pH-dependent equilibria in seawater. These are collated in a dict called `concTotals`, to conveniently pass this information into the solver functions. The dict can be generated using the function of the same name in the `concentrations` module, or you can assemble your own manually.

---

## `concTotals`: dict of concentrations

Assemble a dict of concentrations as required by other Calkulate functions.

**Syntax:**

```python
concTotals = calk.concentrations.concTotals(pSal, totalCarbonate=0,
    totalPhosphate=0, totalSilicate=0, totalAmmonia=0, totalH2Sulfide=0,
    WhichKs=10, WhoseTB=2)
```

The output `concTotals` dict contains the following variables:

*Defined by inputs to the concTotals function:*

  * `C` - `totalCarbonate`.
  * `P` - `totalPhosphate`.
  * `Si` - `totalSilicate`.
  * `NH3` - `totalAmmonia`.
  * `H2S` - `totalH2Sulfide`.

*Estimated from input practical salinity by [PyCO2SYS](https://github.com/mvdh7/PyCO2SYS):*

  * `B` - total borate.
  * `S` - total sulfate.
  * `F` - total fluoride.

The exact functions of salinity used by PyCO2SYS to estimate these concentrations can be set using the inputs `WhichKs`, which is the same as the PyCO2SYS input `K1K2CONSTANTS`, and `WhoseTB`, which corresponds to the PyCO2SYS input `KSO4CONSTANTS`: a value of `1` for `WhoseTB` is the same as options `1` and `2` for `KSO4CONSTANTS`, while a value of `2` is the same as options `3` and `4`.
