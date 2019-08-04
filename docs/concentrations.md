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

Assembles a dict of concentrations as required by other Calkulate functions.

**Syntax:**

```python
concTotals = calk.concentrations.concTotals(pSal, totalCarbonate=0,
    totalPhosphate=0, totalSilicate=0)
```

The output `concTotals` dict contains the following variables:

```python
concTotals = {
# Defined by inputs to the concTotals function:
    'C': totalCarbonate,
    'P': totalPhosphate,
    'Si': totalSilicate,
# Estimated from input practical salinity:
    'B': totalBorate_LKB10(pSal),
    'S': totalSulfate_MR66(pSal),
    'F': totalFluoride_W71(pSal),
}
```

The functions to estimate total concentrations are also contained within the `concentrations` module, along with some alternatives that you could manually implement instead.

---

## Estimates from salinity

The functions available in the `concentrations` module to estimate total concentrations from practical salinity all follow the same syntax:

```python
totalSolute = calk.concentrations.totalSolute_REF(pSal)
```

The different options for `Solute` and `REF` (i.e. [literature reference](../references)) are listed below. The value marked as 'default' is the one used by the `concTotals` function.

### Total borate

`totalBorate` $= [\ce{B(OH)3}] + [\ce{B(OH)4−}]$ in mol/kg-sw.

**Options:**

  * `totalBorate_LKB10`: Lee et al., 2010 [[LKB10](../references/#LKB10)] (default);
  * `totalBorate_U74`: Uppström, 1974 [[U74](../references/#U74)].

### Total fluoride

`totalFluoride` $= [\ce{HF}] + [\ce{F−}]$ in mol/kg-sw.

**Options:**

  * `totalFluoride_W71`: Warner, 1971 [[W71](../references/#W71)] (default).

### Total sulfate

`totalSulfate` $= [\ce{HSO4−}] + [\ce{SO4^2−}]$ in mol/kg-sw.

**Options:**

  * `totalSulfate_MR66`: Morris and Riley, 1966 [[MR66](../references/#MR66)] (default).
