<script type="text/x-mathjax-config">
MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}});
MathJax.Ajax.config.path["mhchem"] =
  "https://cdnjs.cloudflare.com/ajax/libs/mathjax-mhchem/3.3.2";
MathJax.Hub.Config({TeX: {extensions: ["[mhchem]/mhchem.js"]}});
</script><script src='https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-MML-AM_CHTML' async></script>

# Alkalinity solvers

`.solve` contains functions to determine total alkalinity from titration data using several different methods:

  * `.complete` - **recommended**; a least-squares alkalinity and EMF<sub>0</sub> solver using the full expression for total alkalinity defined by Dickson (1981);
  * `.DAA03` - a least-squares alkalinity and $f$ solver following Dickson et al. (2003);
  * `.Dickson1981` - a least-squares alkalinity, dissolved inorganic carbon and $f$ solver following Dickson (1981);
  * `.halfGran` - a modified Gran plot alkalinity and EMF<sub>0</sub> solver following Humphreys (2015).

## Common syntax

**Syntax:**

```python
AT_emf0 = calk.solve.solver(Macid, emf, tempK, Msamp, Cacid, XT, KXF)
```

**Inputs:**

All inputs have [their usual meanings]().

**Output:**

  * `AT_emf0` - full output from the [`scipy.optimize.least_squares`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html) solver:
    * `AT_emf0['x'][0]` - best-fit total alkalinity in mol·kg-sw<sup>−1</sup>;
    * `AT_emf0['x'][1]` - best-fit EMF<sub>0</sub> in mV.

<hr />

## Auxiliary functions

### `.emf2h` - convert EMF to pH

Converts cell potential (EMF) to hydrogen ion concentration using the Nernst equation.

**Syntax:**

```python
h = calk.solve.emf2h(emf, emf0, tempk)
```

**Inputs:**

  * `emf` - cell potential in mV;
  * `emf0` - cell zero potential in mV;
  * `tempk` - seawater temperature in K.

**Output:**

  * `h` - hydrogen ion concentration in mol·kg-sw<sup>−1</sup>.

<br />

### `.h2emf` - convert pH to EMF

Converts hydrogen ion concentration to cell potential (EMF) using the Nernst equation.

**Syntax:**

```python
emf = calk.solve.h2emf(h, emf0, tempk)
```

**Inputs:**

  * `h` - hydrogen ion concentration in mol·kg-sw<sup>−1</sup>;
  * `emf0` - cell zero potential in mV;
  * `tempk` - seawater temperature in K.

**Output:**

  * `emf` - cell potential in mV.

<br />

### `.f2demf0` - $F_1$ to EMF<sub>0</sub>

Converts the $F_1$ function to an estimate of EMF<sub>0</sub>.

**Syntax:**

```python
demf0 = calk.solve.f2demf0(tempK, f)
```

**Inputs:**

  * `tempk` - seawater temperature in K;
  * `f` - $F_1$ function.

**Output:**

  * `demf0` - estimates of EMF<sub>0</sub>.
