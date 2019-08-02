<script type="text/x-mathjax-config">
MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}});
MathJax.Ajax.config.path["mhchem"] =
  "https://cdnjs.cloudflare.com/ajax/libs/mathjax-mhchem/3.3.2";
MathJax.Hub.Config({TeX: {extensions: ["[mhchem]/mhchem.js"]}});
</script><script src='https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-MML-AM_CHTML' async></script>

# Names, dimensions and units

Within Python, the import convention for Calkulate is:

```python
import calkulate as calk
```

The names of variables used and created by Calkulate functions are consistent across the package in terms of both dimensions and units.

General rules for units, unless otherwise specified (and sometimes it is!):

  * Concentrations are in mol per kg of solution (mol/kg-sw for seawater);
  * Electric potentials are in mV;
  * Masses are in kg;
  * pH is on the Free scale i.e. $-log_{10}([\ce{H+}])$
  * Temperatures are in K in Calkulate functions, but °C in input text files;
  * Volumes are in ml.

## Arrays

NumPy arrays of dimension $(n,)$ where $n$ is the number of data points (i.e. acid addition steps) for a given titration:

  * `emf` - electric potential (EMF) measured across the titration cell in mV;
  * `massAcid` - mass of acid added to the titration cell in kg;
  * `pH` - pH at each point in the titration;
  * `tempK` - temperature in the titration cell in K;
  * `volAcid` - volume of acid added to the titration cell in ml.

## Constants

Scalar constants, each of which may have a different value for each titration:

  * `alk` - measured total alkalinity of the original sample in mol/kg-sw;
  * `alkCert` - certified total alkalinity of a reference material in mol/kg-sw;
  * `buretteCorrection` - multiplicative correction factor for the burette volume reported in titration data files (dimensionless);
  * `concAcid` - acid titrant concentration in mol/kg;
  * `CT` - dissolved inorganic carbon in mol/kg-sw;
  * `emf0` - EMF° in mV;
  * `massSample` - mass of sample before any acid addition in kg;
  * `pSal` - practical salinity of the sample (dimensionless);
  * `PT` - phosphate in mol/kg-sw;
  * `SiT` - silicate in mol/kg-sw;
  * `tempKForce` - temperature of the titration in K, if required to override that recorded in a titration data file;
  * `volSample` - volume of sample before any acid addition in ml.

## Dicts

A number of Calkulate functions take two dicts as inputs. These contain:

  1. `XT` - the total concentrations of chemical species that are variable over the titration pH range (for more details see the [concentrations](concentrations) module);
  2. `KXF` - the stoichiometric dissociation constants for each equilibrium in seawater, evaluated at the titration temperature and sample salinity, on the Free pH scale (for more details see the [dissociation](dissociation) module).

## Guesstimates

Calkulate estimates values for quite a few of its variables along the way to calculating the final results. These intermediate guesses have `Guess` appended to the corresponding variable names given in the lists above.
