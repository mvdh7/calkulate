# Solute concentrations

Functions of salinity to estimate the total concentrations of the dissolved components of seawater.


<hr />


## `.XT`: list of concentrations

Assembles a list of concentrations in the order required by other Calkulate functions.

**Syntax:**

```python
XT = calk.conc.XT(psal, CT=0, PT=0, SiT=0)
```

**Inputs:**

  * `psal`: practical salinity;
  * `CT`: dissolved inorganic carbon in mol·kg-sw<sup>−1</sup>;
  * `PT`: total phosphate in mol·kg-sw<sup>−1</sup>;
  * `SiT`: total silicate in mol·kg-sw<sup>−1</sup>.

`CT`, `PT` and `SiT` are optional, and they are assigned values of zero if not declared.

**Output:**

  * `XT`: list of concentrations in the order: `[AT, CT, BT, ST, FT, PT, SiT]` in mol·kg-sw<sup>−1</sup>.

`CT`, `PT` and `SiT` are assigned by the user input, or given values of zero if not specified. `BT`, `ST` and `FT` are estimated from `psal` using the other functions in this module. `AT` is returned with a value of `None`.


<hr />


## `.BT`: total borate

Estimates total borate, i.e. [B(OH)<sub>3</sub>] + [B(OH)<sub>4</sub><sup>−</sup>], from practical salinity following Lee et al. (2010).

**Syntax:**

```python
BT = calk.conc.BT(psal)
```

**Input:**

  * `psal`: practical salinity.

**Output:**

  * `BT`: total borate in mol·kg-sw<sup>−1</sup>.


## `.FT`: total fluoride

Estimates total fluoride, i.e. [HF] + [F<sup>−</sup>], from practical salinity following **?????**.

**Syntax:**

```python
FT = calk.conc.FT(psal)
```

**Input:**

  * `psal`: practical salinity.

**Output:**

  * `FT`: total fluoride in mol·kg-sw<sup>−1</sup>.


## `.ST`: total sulfate

Estimates total sulfate, i.e. [HSO<sub>4</sub><sup>−</sup>] + [SO<sub>4</sub><sup>2−</sup>], from practical salinity following **?????**. Note that the concentration of undissociated H<sub>2</sub>SO<sub>4</sub> is assumed to be zero.

**Syntax:**

```python
ST = calk.conc.ST(psal)
```

**Input:**

  * `psal`: practical salinity.

**Output:**

  * `ST`: total sulfate in mol·kg-sw<sup>−1</sup>.
