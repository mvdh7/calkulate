<script type="text/x-mathjax-config">
MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}});
MathJax.Ajax.config.path["mhchem"] =
  "https://cdnjs.cloudflare.com/ajax/libs/mathjax-mhchem/3.3.2";
MathJax.Hub.Config({TeX: {extensions: ["[mhchem]/mhchem.js"]}});
</script><script src='https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-MML-AM_CHTML' async></script>

# Names, dimensions and units

The import convention for Calkulate is:

```python
import calkulate as calk
```

The names of variables used and created by Calkulate functions are consistent across the package in terms of both dimensions and units.

General rules for units, unless otherwise specified (and beware - sometimes it is!):

  * Concentrations are in mol per kg of solution (mol/kg-sw for seawater);
  * Electric potentials are in mV;
  * Masses are in kg;
  * pH and equilibrium constants are on the Free scale, i.e. $\text{pH} = -\log_{10}([\ce{H+}])$;
  * Temperatures are in K in Calkulate functions, but °C in input text files;
  * Volumes are in ml.

## Arrays

NumPy arrays of dimension $(n,)$ where $n$ is the number of data points (i.e. acid addition steps) for a given titration:

  * `emf` - electric potential (EMF) measured across the titration cell solution in mV;
  * `massAcid` - mass of acid added to the sample in kg;
  * `pH` - solution pH at each point in the titration;
  * `tempK` - temperature of the solution in the titration cell in °C;
  * `tempK` - temperature of the solution in the titration cell in K;
  * `volAcid` - volume of acid added to the sample in ml.

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
  * `totalBorate` - total borate in the sample in mol/kg:
    - $[\ce{B(OH)3}] + [\ce{B(OH)4-}]$
  * `totalCarbonate` - total carbonate (dissolved inorganic carbon) in the sample in mol/kg:
    - $[\ce{CO2(aq)}] + [\ce{HCO3-}] + [\ce{CO3^2-}]$;
  * `totalFluoride` - total fluoride in the sample in mol/kg:
    - $[\ce{HF}] + [\ce{F-}]$;
  * `totalPhosphate` - total phosphate in the sample in mol/kg:
    - $[\ce{H3PO4}] + [\ce{H2PO4-}] + [\ce{HPO4^2-}] + [\ce{PO4^3-}]$;
  * `totalSilicate` - total silicate in the sample in mol/kg:
    - $[\ce{Si(OH)4}] + [\ce{SiO(OH)3-}]$
  * `totalSulfate` - total sulfate in the sample in mol/kg:
    - $[\ce{HSO4-}] + [\ce{SO4^2-}]$
  * `volSample` - volume of sample before any acid addition in ml.

## Dicts

A number of Calkulate functions take two dicts as inputs. These are, and contain:

  1. `concTotals` - the constant total concentrations of chemical species that are individually variable over the titration pH range (for more details see the [total concentrations](../concentrations) documentation);
  2. `eqConstants` - the stoichiometric dissociation constants for each equilibrium in seawater, evaluated at the titration temperature and sample salinity, on the Free pH scale (for more details see the [dissociation constants](../dissociation) documentation).

## Guesstimates

Calkulate estimates values for quite a few of its variables along the way to calculating the final results (e.g. to provide initial estimates for the least-squares solvers). These intermediate guesses have `Guess` appended to the corresponding variable names given in the lists above.
