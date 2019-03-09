# Alkalinity solvers

Contains functions to determine total alkalinity from titration data using several different methods.


<hr />

## .EMF2H - convert EMF to pH

Converts cell potential (EMF) to hydrogen ion concentration using the Nernst equation.

**Syntax:**

```python
H = calk.solve.EMF2H(EMF, EMF0, tempK)
```

**Inputs:**

  * `EMF`: cell potential in mV;
  * `EMF0`: cell zero potential in mV;
  * `tempK`: seawater temperature in K.

**Output:**

  * `H`: hydrogen ion concentration in mol·kg-sw<sup>−1</sup>.


## .H2EMF - convert pH to EMF

Converts hydrogen ion concentration to cell potential (EMF) using the Nernst equation.

**Syntax:**

```python
EMF = calk.solve.H2EMF(H, EMF0, tempK)
```

**Inputs:**

  * `H`: hydrogen ion concentration in mol·kg-sw<sup>−1</sup>;
  * `EMF0`: cell zero potential in mV;
  * `tempK`: seawater temperature in K.

**Output:**

  * `EMF`: cell potential in mV.
