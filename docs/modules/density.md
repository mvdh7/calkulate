# Fluid densities

`.density` contains functions of temperature and salinity that estimate the densities of [1] seawater and [2] mixtures of hydrochloric acid (HCl) and sodium chloride (NaCl).


<hr />

## `.sw` - seawater at 1 atm

Calculates the density of seawater as a function of its temperature and practical salinity at a pressure of 1 atm, following Millero and Poisson (1981).

**Syntax:**

```python
rho_sw = calk.density.sw(tempK, psal)
```

**Inputs:**

  * `tempK` - seawater temperature in K;
  * `psal` - practical salinity.

**Output:**

  * `rho_sw` - seawater density in kg·dm<sup>−3</sup>.


<hr />

## `.acid` - 0.1 M HCl in 0.6 M NaCl

Calculates the density of a mixed solution of 0.1 M HCl and 0.6 M NaCl as a function of its temperature, at a pressure of 1 atm.

The function uses a second-order polynomial fit through a series of temperature-density pairs, from 288.15 to 308.15 K, calculated using [E-AIM](http://www.aim.env.uea.ac.uk/aim/density/density_electrolyte.php) (Clegg and Wexler, 2011a, 2011b) with option 'rho, at the total solute mass fraction' and the concentrations:

  * [H<sup>+</sup>] = 0.1 mol·dm<sup>−3</sup>;
  * [Na<sup>+</sup>] = 0.6 mol·dm<sup>−3</sup>;
  * [Cl<sup>−</sup>] = 0.7 mol·dm<sup>−3</sup>.

This represents a 0.1 mol·dm<sup>−3</sup> HCl titrant mixed with NaCl, with total ionic strength like that of typical open ocean seawater (Dickson et al., 2003).

**Syntax:**

```python
rho_acid = calk.density.acid(tempK)
```

**Input:**

  * `tempK` - acid temperature in K.

**Output:**

  * `rho_acid` - acid density in kg·dm<sup>−3</sup>.

<hr />

## `.acid25` - HCl in NaCl at 25 °C

Calculates the density of a mixed solution of HCl and NaCl as a function of its composition, at a pressure of 1 atm and temperature of 25 °C, following Dickson et al. (2007).

Note that the result from this function does not quite agree with the check value provided by Dickson et al. (2007). They declare it should give 1.02056 kg·dm<sup>−3</sup> for HCl and NaCl concentrations of 0.2 and 0.5 mol·kg-H<sub>2</sub>O<sup>−1</sup> respectively. However, this function returns 1.02035 kg·dm<sup>−3</sup> instead. The reason for this discrepancy is unclear.

**Syntax:**

```python
rho25 = calk.density.acid25(mHCl, mNaCl)
```

**Input:**

  * `mHCl` - HCl concentration in mol·kg-H<sub>2</sub>O<sup>−1</sup>;
  * `mNaCl` - NaCl concentration in mol·kg-H<sub>2</sub>O<sup>−1</sup>.

**Output:**

  * `rho25` - acid density in kg·dm<sup>−3</sup>.
