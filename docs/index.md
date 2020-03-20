# <img src="img/logo_transparent.png" style="vertical-align:sub" width="105px" /> v2.2.0

Calkulate is a Python 3.6+ package for determining total alkalinity from seawater titration data.

## Installation

If you're using Conda, first create a new environment with Python 3.6+, NumPy 1.15+, and SciPy 1.1+ - or you can just allow Pip to install these dependencies for you. Other similar versions are probably fine, but untested. Then, to install:

    pip install calkulate

To upgrade in existing installation when there is a new release:

    pip install calkulate --upgrade --no-cache-dir

Since v2.2.0, Calkulate also requires [PyCO2SYS](https://github.com/mvdh7/PyCO2SYS) v1.1.1 or greater. This will be automatically installed by Pip if you don't have it.

## Get Calkulating!

This documentation is intended to present a broad overview of how Calkulate works, rather than provide a detailed syntactic reference for every constituent module and function.

Within Python, the import convention is:

```python
import calkulate as calk
```

To quickly get started with some analysis, jump straight to the [workflow examples](../examples/compare-all-solvers). These examples illustrate the different parts of Calkulate that you may need to use, and provide a framework that you can quickly modify to suit your own requirements.

To find out more about the principles by which Calkulate calculates things, take a look at the other parts of the documentation. These explain:

  * [Variables and conventions](conventions): the conventions for naming and defining the input and output variables;
  * [Data import](io): how to import titration data to work with;
  * [Total concentrations](concentrations): estimating total concentrations from salinity, or defining your own;
  * [Equilibrium constants](dissociation): defining stoichiometric equilibrium constants for the alkalinity solvers;
  * [Alkalinity solvers](solvers): the different methods by which alkalinity can be estimated from titration data;
  * [Titrant calibration](calibration): how to calibrate the acid titrant's concentration;
  * [Version history](versions): changes from version to version; and
  * [Literature references](references): a key to the codes used for citations.

## Citation

A paper describing Calkulate v2 is in preparation. For now, if you use any version of Calkulate in your research, please cite it as:

> Humphreys, M. P. (2015). "Calculating seawater total alkalinity from open-cell titration data using a modified Gran plot technique," in *Measurements and Concepts in Marine Carbonate Chemistry* (PhD Thesis, Ocean and Earth Science, University of Southampton, UK), 25–44.

But please do check back here for any updates first (or [get in touch](https://mvdh.xyz/contact/))!

## Acknowledgements

Calkulate is maintained by [Dr Matthew P. Humphreys](https://mvdh.xyz) at the NIOZ Royal Netherlands Institute for Sea Research (Texel, the Netherlands), and Ruth Matthews at the University of East Anglia (Norwich, UK).

Its ongoing development has been indirectly funded by the [Natural Environment Research Council](https://nerc.ukri.org/) (NERC, UK) and the [Dutch Research Council](https://www.nwo.nl/en) (NWO, the Netherlands).

<!--
  * Originally through a PhD studentship to Matthew P. Humphreys (NE/J500112/1),
  * Followed by *CaNDyFloSS: Carbon and Nutrient Dynamics and Fluxes over Shelf Systems* (NE/K00185X/1) and *RAGNARoCC: Radiatively active gases from the North Atlantic Region and Climate Change* (NE/K002546/1),
  * Then through *NSFGEO-NERC: A Thermodynamic Chemical Speciation Model for the Oceans, Seas, and Estuaries* (NE/P012361/1).
-->

## License

The entirety of Calkulate is licensed under the [GNU General Public License version 3 (GPLv3)](https://www.gnu.org/licenses/gpl-3.0.en.html).
