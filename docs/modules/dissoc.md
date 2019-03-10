# Dissociation constants

Functions of temperature and salinity that estimate the stoichiometric dissociation constants required to model seawater equilibria.


<!-- Add chemistry with https://mhchem.github.io/MathJax-mhchem/ -->

<hr />


## .KX_F - list of constants

Assembles a list of dissociation constants on the Free pH scale in the order required by other Calkulate functions.

**Syntax:**

```python
KX = calk.dissoc.KX_F(tempK, psal, ST, FT)
```

**Inputs:**

  * `tempK`: titre temperature in K;
  * `psal`: practical salinity;
  * `ST`: total sulfate concentration in mol·kg-sw<sup>−1</sup>;
  * `FT`: total fluoride concentration in mol·kg-sw<sup>−1</sup>.

**Output:**

  * `KX`: list of dissociation constants in the order: `[KC1, KC2, KB, KH2O, KHSO4, KHF, KP1, KP2, KP3, KSi]`, all on the Free pH scale.

The dissociation constants are evaluated using the other functions in this module.


<hr />

## .Istr - ionic strength

Estimates ionic strength from salinity following **?????**.

**Syntax:**

```python
Istr = calk.dissoc.Istr(psal)
```

**Input:**

  * `psal`: practical salinity.

**Output:**

  * `Istr`: ionic strength in mol·kg-sw<sup>−1</sup>.


<hr />

## .KC_T_LDK00 - carbonic acid

Estimates the carbonic acid stoichiometric dissociation constants K1 and K2 from temperature and salinity on the Total pH scale following Lueker et al. (2000). Valid for temperature from 2 to 35 degC and practical salinity from 19 to 43.

**Syntax:**

```python
K1, K2 = calk.dissoc.KC_T_LDK00(tempK, psal)
```

**Inputs:**

  * `tempK`: temperature in K;
  * `psal`: practical salinity.

**Outputs:**

  * `K1`: first dissociation constant for carbonic acid on the Total pH scale in mol·kg-sw<sup>−1</sup>;
  * `K2`: second dissociation constant for carbonic acid on the Total pH scale in mol·kg-sw<sup>−1</sup>.


## .KB_T_D90a - boric acid

Estimates the boric acid stoichiometric dissociation constant KB from temperature and salinity on the Total pH scale following Dickson (1990a). Valid for temperature from 0 to 45 degC and practical salinity from 5 to 45.

**Syntax:**

```python
KB = calk.dissoc.KB_T_D90a(tempK, psal)
```

**Inputs:**

  * `tempK`: temperature in K;
  * `psal`: practical salinity.

**Outputs:**

  * `KB`: dissociation constant for boric acid on the Total pH scale in mol·kg-sw<sup>−1</sup>.


## .KH2O_T_DSC07 - water

Estimates the water stoichiometric dissociation constant KH2O from temperature and salinity on the Total pH scale following Dickson et al. (2007).

**Syntax:**

```python
KH2O = calk.dissoc.KH2O_T_DSC07(tempK, psal)
```

**Inputs:**

  * `tempK`: temperature in K;
  * `psal`: practical salinity.

**Outputs:**

  * `KH2O`: dissociation constant for water acid on the Total pH scale in mol·kg-sw<sup>−1</sup>.
