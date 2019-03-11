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

Estimates the carbonic acid stoichiometric dissociation constants K1 and K2 from temperature and salinity on the Total pH scale following Lueker et al. (2000). Valid for temperature from 2 to 35 °C and practical salinity from 19 to 43.

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

Estimates the boric acid stoichiometric dissociation constant KB from temperature and salinity on the Total pH scale following Dickson (1990a). Valid for temperature from 0 to 45 °C and practical salinity from 5 to 45.

**Syntax:**

```python
KB = calk.dissoc.KB_T_D90a(tempK, psal)
```

**Inputs:**

  * `tempK`: temperature in K;
  * `psal`: practical salinity.

**Output:**

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

**Output:**

  * `KH2O`: dissociation constant for water acid on the Total pH scale in mol·kg-sw<sup>−1</sup>.


## .KHSO4_F_D90b - bisulfate

Estimates the bisulfate stoichiometric dissociation constant KHSO4 from temperature and salinity on the Free pH scale following Dickson (1990b). Valid for temperature from 0 to 45 °C and practical salinity from 5 to 45.

**Syntax:**

```python
KHSO4 = calk.dissoc.KHSO4_F_D90b(tempK, psal)
```

**Inputs:**

  * `tempK`: temperature in K;
  * `psal`: practical salinity.

**Output:**

  * `KHSO4`: dissociation constant for bisulfate on the Free pH scale in mol·kg-sw<sup>−1</sup>.


## .KHF_T_PF87 - hydrogen fluoride

Estimates the hydrogen fluoride stoichiometric dissociation constant KHF from temperature and salinity on the Total pH scale following Perez and Fraga (1987). Valid for temperature from 9 to 33 °C and practical salinity from 10 to 40.

**Syntax:**

```python
KHF = calk.dissoc.KHF_T_PF87(tempK, psal)
```

**Inputs:**

  * `tempK`: temperature in K;
  * `psal`: practical salinity.

**Output:**

  * `KHF`: dissociation constant for hydrogen fluoride on the Total pH scale in mol·kg-sw<sup>−1</sup>.


## .KHF_F_DR79 - hydrogen fluoride

Estimates the hydrogen fluoride stoichiometric dissociation constant KHF from temperature and salinity on the Free pH scale following Dickson and Riley (1979). Valid for temperature from 5 to 35 °C and practical salinity from 10 to 48.

**Syntax:**

```python
KHF = calk.dissoc.KHF_F_DR79(tempK, psal)
```

**Inputs:**

  * `tempK`: temperature in K;
  * `psal`: practical salinity.

**Output:**

  * `KHF`: dissociation constant for hydrogen fluoride on the Free pH scale in mol·kg-sw<sup>−1</sup>.


## .KP_T_DSC07 - phosphoric acid

Estimates the phosphoric acid stoichiometric dissociation constants KP1, KP2 and KP3 from temperature and salinity on the Total pH scale following Dickson et al. (2007).

**Syntax:**

```python
KP1, KP2, KP3 = calk.dissoc.KP_T_DSC07(tempK, psal)
```

**Inputs:**

  * `tempK`: temperature in K;
  * `psal`: practical salinity.

**Outputs:**

  * `KP1`: first dissociation constant for phosphoric acid on the Total pH scale in mol·kg-sw<sup>−1</sup>;
  * `KP2`: second dissociation constant for phosphoric acid on the Total pH scale in mol·kg-sw<sup>−1</sup>;
  * `KP3`: third dissociation constant for phosphoric acid on the Total pH scale in mol·kg-sw<sup>−1</sup>.


## .KSi_T_M95 - silicic acid

Estimates the silicic acid stoichiometric dissociation constant from temperature and salinity on the Total pH scale following Millero (1995).

**Syntax:**

```python
KSi = calk.dissoc.KSi_T_M95(tempK, psal)
```

**Inputs:**

  * `tempK`: temperature in K;
  * `psal`: practical salinity.

**Outputs:**

  * `KSi`: dissociation constant for silicic acid on the Total pH scale in mol·kg-sw<sup>−1</sup>.


## .KNH4_X_BJJL08 - ammonium

Estimates the ammonium stoichiometric dissociation constant from temperature and salinity on **an unknown** pH scale following Bell et al. (2008).

**Syntax:**

```python
KNH4 = calk.dissoc.KNH4_X_BJJL08(tempK, psal)
```

**Inputs:**

  * `tempK`: temperature in K;
  * `psal`: practical salinity.

**Outputs:**

  * `KNH4`: dissociation constant for ammonium on **an unknown** pH scale in mol·kg-sw<sup>−1</sup>.


## .K2AMP_S_BE86 - 2-aminopyridine

Estimates the 2-aminopyridine stoichiometric dissociation constant from temperature and salinity on the seawater pH scale following Bates and Erickson (1986). Valid for temperature from 5 to 40 °C and practical salinity from 30 to 40.

**Syntax:**

```python
K2AMP = calk.dissoc.K2AMP_S_BE86(tempK, psal)
```

**Inputs:**

  * `tempK`: temperature in K;
  * `psal`: practical salinity.

**Outputs:**

  * `K2AMP`: dissociation constant for 2-aminopyridine on the seawater pH scale in mol·kg-sw<sup>−1</sup>.
