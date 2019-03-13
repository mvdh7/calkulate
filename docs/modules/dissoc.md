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


# Dissociation constants

Functions of temperature and salinity that estimate the stoichiometric dissociation constants required to model seawater equilibria.


<!-- Add chemistry with https://mhchem.github.io/MathJax-mhchem/ -->

<hr />


## `.KX_F`: list of constants

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

## `.Istr`: ionic strength

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

# Carbonic acid

After dissolving into seawater, $\ce{CO2}$ reacts to form carbonic acid ($\ce{H2CO3}$), which totally dissociates into bicarbonate and carbonate ions ($\ce{HCO3-}$ and $\ce{CO3^2-}$). The relevant equilibria are:

$$\ce{CO2(aq) + H2O <=>[$K_1$] HCO3- + H+}$$

$$\ce{HCO3- <=>[$K_2$] CO3^2- + H+}$$

with the stoichiometric dissociation constants:

$$K_1^\* = \frac{[\ce{HCO3-}] [\ce{H+}]}{[\ce{CO2(aq)}]}$$

$$K_2^\* = \frac{[\ce{CO3^2-}] [\ce{H+}]}{[\ce{HCO3-}]}$$

## `.KC_T_LDK00`: Lueker et al. (2000)

Estimates the carbonic acid stoichiometric dissociation constants $K_1^\*$ and $K_2^\*$ from temperature and salinity on the Total pH scale following Lueker et al. (2000).

**Validity:**

  * Temperature: from 2 to 35 °C;
  * Practical salinity: from 19 to 43.

**Syntax:**

```python
K1, K2 = calk.dissoc.KC_T_LDK00(tempK, psal)
```

**Inputs:**

  * `tempK`: temperature in K;
  * `psal`: practical salinity.

**Outputs:**

  * `K1`: first stoichiometric dissociation constant for carbonic acid ($K_1^\*$), on the Total pH scale, in mol·kg-sw<sup>−1</sup>;
  * `K2`: second stoichiometric dissociation constant for carbonic acid ($K_2^\*$), on the Total pH scale, in mol·kg-sw<sup>−1</sup>.

<hr />

# Boric acid

The equilibrium reaction is:

$$\ce{B(OH)3 + H2O <=>[$K_\text{B}$] B(OH)4- + H+}$$

with the stoichiometric dissociation constant:

$$K_\text{B}^\* = \frac{[\ce{B(OH)4-}] [\ce{H+}]}{[\ce{B(OH)3}]}$$

## `.KB_T_D90a`: Dickson (1990a)

Estimates the boric acid stoichiometric dissociation constant $K_\text{B}^\*$ from temperature and salinity on the Total pH scale following Dickson (1990a).

**Validity:**

  * Temperature: from 0 to 45 °C;
  * Practical salinity: from 5 to 45.

**Syntax:**

```python
KB = calk.dissoc.KB_T_D90a(tempK, psal)
```

**Inputs:**

  * `tempK`: temperature in K;
  * `psal`: practical salinity.

**Output:**

  * `KB`: stoichiometric dissociation constant for boric acid ($K_\text{B}^\*$), on the Total pH scale, in mol·kg-sw<sup>−1</sup>.

<hr />

# Water

The equilibrium reaction is:

$$\ce{H2O <=>[$K_\text{w}$] OH- + H+}$$

with the stoichiometric dissociation constant:

$$K_\text{w}^\* = [\ce{OH-}] [\ce{H+}]$$

## `.KH2O_T_DSC07`: Dickson et al. (2007)

Estimates the water stoichiometric dissociation constant $K_\text{w}^\*$ from temperature and salinity on the Total pH scale following Dickson et al. (2007).

**Syntax:**

```python
KH2O = calk.dissoc.KH2O_T_DSC07(tempK, psal)
```

**Inputs:**

  * `tempK`: temperature in K;
  * `psal`: practical salinity.

**Output:**

  * `KH2O`: stoichiometric dissociation constant for water ($K_\text{w}^\*$), on the Total pH scale, in (mol·kg-sw<sup>−1</sup>)<sup>2</sup>.

<hr />

# Bisulfate

The equilibrium reaction is:

$$\ce{HSO4- <=>[$K_\text{S}$] SO4^2- + H+}$$

with the stoichiometric dissociation constant:

$$K_\text{S}^\* = \frac{[\ce{SO4^2-}] [\ce{H+}]}{[\ce{HSO4-}]}$$

## `.KHSO4_F_D90b`: Dickson (1990b)

Estimates the bisulfate stoichiometric dissociation constant $K_\text{w}^\*$ from temperature and salinity on the Free pH scale following Dickson (1990b).

**Validity:**

  * Temperature: from 0 to 45 °C;
  * Practical salinity: from 5 to 45.

**Syntax:**

```python
KHSO4 = calk.dissoc.KHSO4_F_D90b(tempK, psal)
```

**Inputs:**

  * `tempK`: temperature in K;
  * `psal`: practical salinity.

**Output:**

  * `KHSO4`: stoichiometric dissociation constant for bisulfate ($K_\text{w}^\*$), on the Free pH scale, in mol·kg-sw<sup>−1</sup>.

<hr />

# Hydrofluoric acid

The equilibrium reaction is:

$$\ce{HF <=>[$K_\text{HF}$] F- + H+}$$

with the stoichiometric dissociation constant:

$$K_\text{HF}^\* = \frac{[\ce{F-}] [\ce{H+}]}{[\ce{HF}]}$$

## `.KHF_T_PF87`: Perez and Fraga (1987)

Estimates the hydrogen fluoride stoichiometric dissociation constant $K_\text{HF}^\*$ from temperature and salinity on the Total pH scale following Perez and Fraga (1987).

**Validity:**

  * Temperature: from 9 to 33 °C;
  * Practical salinity: from 10 to 40.

**Syntax:**

```python
KHF = calk.dissoc.KHF_T_PF87(tempK, psal)
```

**Inputs:**

  * `tempK`: temperature in K;
  * `psal`: practical salinity.

**Output:**

  * `KHF`: stoichiometric dissociation constant for hydrogen fluoride ($K_\text{HF}^\*$), on the Total pH scale, in mol·kg-sw<sup>−1</sup>.


## `.KHF_F_DR79`: Dickson and Riley (1979)

Estimates the hydrogen fluoride stoichiometric dissociation constant $K_\text{HF}^\*$ from temperature and salinity on the Free pH scale following Dickson and Riley (1979).

**Validity:**

  * Temperature: from 5 to 35 °C;
  * Practical salinity: from 10 to 48.

**Syntax:**

```python
KHF = calk.dissoc.KHF_F_DR79(tempK, psal)
```

**Inputs:**

  * `tempK`: temperature in K;
  * `psal`: practical salinity.

**Output:**

  * `KHF`: stoichiometric dissociation constant for hydrogen fluoride ($K_\text{HF}^\*$), on the Free pH scale, in mol·kg-sw<sup>−1</sup>.

<hr />

# Phosphoric acid

The equilibrium reactions are:

$$\ce{H3PO4 <=>[$K_\text{P1}$] H2PO4- + H+}$$

$$\ce{H2PO4- <=>[$K_\text{P2}$] HPO4^2- + H+}$$

$$\ce{HPO4^2- <=>[$K_\text{P3}$] PO4^3- + H+}$$

with the stoichiometric dissociation constants:

$$K_\text{P1}^\* = \frac{[\ce{H2PO4-}] [\ce{H+}]}{[\ce{H3PO4}]}$$

$$K_\text{P2}^\* = \frac{[\ce{HPO4^2-}] [\ce{H+}]}{[\ce{H2PO4-}]}$$

$$K_\text{P3}^\* = \frac{[\ce{PO4^3-}] [\ce{H+}]}{[\ce{HPO4^2-}]}$$

## `.KP_T_DSC07 `: Dickson et al. (2007)

Estimates the phosphoric acid stoichiometric dissociation constants $K_\text{P1}^\*$, $K_\text{P2}^\*$ and $K_\text{P3}^\*$ from temperature and salinity on the Total pH scale following Dickson et al. (2007).

**Syntax:**

```python
KP1, KP2, KP3 = calk.dissoc.KP_T_DSC07(tempK, psal)
```

**Inputs:**

  * `tempK`: temperature in K;
  * `psal`: practical salinity.

**Outputs:**

  * `KP1`: first stoichiometric dissociation constant for phosphoric acid ($K_\text{P1}^\*$), on the Total pH scale, in mol·kg-sw<sup>−1</sup>;
  * `KP2`: second stoichiometric dissociation constant for phosphoric acid ($K_\text{P2}^\*$), on the Total pH scale, in mol·kg-sw<sup>−1</sup>;
  * `KP3`: third stoichiometric dissociation constant for phosphoric acid ($K_\text{P3}^\*$), on the Total pH scale, in mol·kg-sw<sup>−1</sup>.

<hr />

# Orthosilicic acid

The equilibrium reaction is:

$$\ce{Si(OH)4 <=>[$K_\text{Si}$] Si(OH)3O- + H+}$$

with the stoichiometric dissociation constant:

$$K_\text{Si}^\* = \frac{[\ce{Si(OH)3O-}] [\ce{H+}]}{[\ce{Si(OH)4}]}$$

## `.KSi_T_M95`: Millero (1995)

Estimates the orthosilicic acid stoichiometric dissociation constant $K_\text{Si}^\*$ from temperature and salinity on the Total pH scale following Millero (1995).

**Syntax:**

```python
KSi = calk.dissoc.KSi_T_M95(tempK, psal)
```

**Inputs:**

  * `tempK`: temperature in K;
  * `psal`: practical salinity.

**Outputs:**

  * `KSi`: stoichiometric dissociation constant for silicic acid ($K_\text{Si}^\*$), on the Total pH scale, in mol·kg-sw<sup>−1</sup>.

<hr />

<!--

# Ammonium

The equilibrium reaction is:

$$\ce{NH3 <=>[$K_\text{Si}$] Si(OH)3O- + H+}$$

with the stoichiometric dissociation constant:

$$K_\text{Si}^\* = \frac{[\ce{Si(OH)3O-}] [\ce{H+}]}{[\ce{Si(OH)4}]}$$

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

-->
