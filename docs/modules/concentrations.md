<script type="text/x-mathjax-config">
MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}});
MathJax.Ajax.config.path["mhchem"] =
  "https://cdnjs.cloudflare.com/ajax/libs/mathjax-mhchem/3.3.2";
MathJax.Hub.Config({
  TeX: {
    extensions: ["[mhchem]/mhchem.js"]
  }
});
</script>
<script src='https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-MML-AM_CHTML' async></script>


# Solute concentrations from salinity

`.concentrations` contains functions of salinity to estimate the total concentrations of the dissolved components of seawater.


<hr />


## `.XT` - dict of concentrations

Assembles a dict of concentrations as required by other Calkulate functions.

**Syntax:**

```python
XT = calk.concentrations.XT(psal, CT=0, PT=0, SiT=0)
```

**Inputs:**

  * `psal` - practical salinity;
  * `CT` - dissolved inorganic carbon in mol·kg-sw<sup>−1</sup>;
  * `PT` - total phosphate in mol·kg-sw<sup>−1</sup>;
  * `SiT` - total silicate in mol·kg-sw<sup>−1</sup>.

`CT`, `PT` and `SiT` are optional, and they are assigned values of zero if not declared.

**Output:**

  * `XT` - dict of concentrations in mol·kg-sw<sup>−1</sup>.

Output fields `C`, `P` and `Si` are assigned by the respective user inputs, or given values of zero if not specified. `B`, `S` and `F` are estimated from `psal` using the other functions in this module.


<hr />


## `.BT` - total borate

Estimates total borate, i.e. $[\ce{B(OH)3}] + [\ce{B(OH)4−}]$, from practical salinity following Lee et al. (2010).

**Syntax:**

```python
BT = calk.concentrations.BT(psal)
```

**Input:**

  * `psal` - practical salinity.

**Output:**

  * `BT` - total borate in mol·kg-sw<sup>−1</sup>.


## `.FT` - total fluoride

Estimates total fluoride, i.e. $[\ce{HF}] + [\ce{F−}]$, from practical salinity following **?????**.

**Syntax:**

```python
FT = calk.concentrations.FT(psal)
```

**Input:**

  * `psal` - practical salinity.

**Output:**

  * `FT` - total fluoride in mol·kg-sw<sup>−1</sup>.


## `.ST` - total sulfate

Estimates total sulfate, i.e. $[\ce{HSO4−}] + [\ce{SO4^2−}]$, from practical salinity following **?????**. The concentration of undissociated $\ce{H2SO4}$ is assumed to be zero.

**Syntax:**

```python
ST = calk.concentrations.ST(psal)
```

**Input:**

  * `psal` - practical salinity.

**Output:**

  * `ST` - total sulfate in mol·kg-sw<sup>−1</sup>.
