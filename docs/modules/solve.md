# Alkalinity solvers

`.solve` contains functions to determine total alkalinity from titration data using several different methods.


<hr />

## `.emf2h` - convert EMF to pH

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


## `.h2emf` - convert pH to EMF

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
