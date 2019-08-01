# Variable naming conventions

The names of variables used and created by Calkulate functions are intended to be consistent across the package.

## Arrays

The following variables are NumPy arrays of dimension *n*, where *n* is the number of data points (i.e. acid addition steps) for a given titration:

  * `emf` - electric potential measured across the titration cell in mV;
  * `massAcid` - mass of acid added to the titration cell in kg;
  * `pH` - Free scale pH at each point in the titration;
  * `tempK` - temperature in the titration cell in K;
  * `volAcid` - volume of acid added to the titration cell in ml.

## Constants

The following variables are scalar constants, which may have different values for each titration:

  * `alk` - measured total alkalinity of the original sample in mol/kg-sw;
  * `alkCert` - certified total alkalinity of a reference material in mol/kg-sw;
  * `buretteCorrection` - multiplicative correction factor for the burette volume reported in titration data files (dimensionless);
  * `concAcid` - acid titrant concentration in mol/kg;
  * `CT` - dissolved inorganic carbon in mol/kg-sw;
  * `emf0` - EMF0 in mV;
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
