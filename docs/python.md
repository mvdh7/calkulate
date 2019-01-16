# calkulate.dens

Functions to estimate the densities of various fluids.



## calkulate.dens.sw

### Description

Calculates the density of seawater as a function of its temperature and practical salinity, at a pressure of 1 atm, following Millero and Poisson (1981).

### Syntax

```python
dens_sw = calk.dens.sw(Tk,S)
```

### Inputs

<table><tr>

<td><strong>Variable</strong></td>
<td><strong>Description</strong></td>
<td><strong>Unit</strong></td>

</tr><tr>

<td><code>Tk</code></td>
<td>Seawater temperature</td>
<td>K</td>

</tr><tr>

<td><code>S</code></td>
<td>Practical salinity</td>
<td>none</td>

</tr></table>

### Output

<table><tr>

<td><strong>Variable</strong></td>
<td><strong>Description</strong></td>
<td><strong>Unit</strong></td>

</tr><tr>

<td><code>dens_sw</code></td>
<td>Seawater density at atmospheric pressure</td>
<td>kg·ℓ<sup>−1</sup></td>

</tr></table>



## calkulate.dens.acid

### Description

Calculates the density of a mixed solution of HCl and NaCl as a function of its temperature, at a pressure of 1 atm.

Uses a second-order polynomial fit through a series of temperature-density pairs, from 288.15 to 308.15 K, calculated using [E-AIM](http://www.aim.env.uea.ac.uk/aim/density/density_electrolyte.php) with option 'rho, at the total solute mass fraction' and the concentrations: [H<sup>+</sup>] = 0.1 M; [Na<sup>+</sup>] = 0.6 mol·ℓ<sup>−1</sup>; [Cl<sup>−</sup>] = 0.7 mol·ℓ<sup>−1</sup>. This represents a 0.1 mol·ℓ<sup>−1</sup> HCl titrant mixed with NaCl, with total ionic strength equal to that of typical open ocean seawater (Dickson et al., 2003).

### Syntax

```python
dens_acid = calk.dens.acid(Tk)
```

### Input

<table><tr>

<td><strong>Variable</strong></td>
<td><strong>Description</strong></td>
<td><strong>Unit</strong></td>

</tr><tr>

<td><code>Tk</code></td>
<td>Acid temperature</td>
<td>K</td>

</tr></table>

### Output

<table><tr>

<td><strong>Variable</strong></td>
<td><strong>Description</strong></td>
<td><strong>Unit</strong></td>

</tr><tr>

<td><code>dens_acid</code></td>
<td>Acid density</td>
<td>kg·ℓ<sup>−1</sup></td>

</tr></table>




# References

Dickson, A. G., Afghan, J. D., and Anderson, G. C. (2003). Reference materials for oceanic CO<sub>2</sub> analysis: a method for the certification of total alkalinity. *Mar. Chem.* 80, 185–197. <a href="https://doi.org/10.1016/S0304-4203(02)00133-0">doi:10.1016/S0304-4203(02)00133-0</a>.

Millero, F. J., and Poisson, A. (1981). International one-atmosphere equation of state of seawater. *Deep-Sea Res. Pt. A* 28, 625–629. <a href="https://doi.org/10.1016/0198-0149(81)90122-9">doi:10.1016/0198-0149(81)90122-9</a>.
