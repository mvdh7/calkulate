# Solute concentrations

Functions of salinity to estimate the total concentrations of the dissolved components of seawater.


## .XT - create concentration list

Assembles a list of concentrations in the order required by other Calkulate functions.

**Syntax:**

```python
XT = calk.conc.XT(S, CT=0, PT=0, SiT=0)
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


## .BT - total borate

Estimates total borate from practical salinity following Lee et al. (2010).
