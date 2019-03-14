<!--<script src='https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-MML-AM_CHTML' async></script>-->

<div style="text-align: center; padding-bottom: 4%">
<img src="img/logo_transparent.png" width="50%" />
</div>

**Calkulate** is a Python 3.6+ package for determining total alkalinity from seawater titration data.


# Installation

If using conda, first create and activate a new environment with Python 3.6+, numpy 1.15+, and scipy 1.1+ (or you can allow pip to install these dependencies for you). Then:

```shell
pip install calkulate
```

Other versions are probably fine, but untested. We use Python 3.6 rather than 3.7 to enable integration with MATLAB.

See the [quick-start guide](quick-start) for more detailed instructions and examples.


# Development status

Calkulate v2.0 is in beta. Tests of the accuracy of its coefficients and equations are underway, so results may change. API may change and functions may be added or removed. Use at your own peril!


# Acknowledgements

Calkulate is maintained by [Dr Matthew P. Humphreys](https://mvdh.xyz) at the Centre for Ocean and Atmospheric Sciences, School of Environmental Sciences, University of East Anglia, Norwich, UK.

Its ongoing development has been funded by the [Natural Environment Research Council](https://nerc.ukri.org/) (NERC, UK):

  * Originally through a PhD studentship to Matthew P. Humphreys (NE/J500112/1),
  * Followed by *CaNDyFloSS: Carbon and Nutrient Dynamics and Fluxes over Shelf Systems* (NE/K00185X/1) and *RAGNARoCC: Radiatively active gases from the North Atlantic Region and Climate Change* (NE/K002546/1),
  * And now through *NSFGEO-NERC: A Thermodynamic Chemical Speciation Model for the Oceans, Seas, and Estuaries* (NE/P012361/1).


# License

<img src="img/1920px-GPLv3_Logo.svg.png" width="25%" />

The entirety of Calkulate is licensed under the [GNU General Public License version 3 (GPLv3)](https://www.gnu.org/licenses/gpl-3.0.en.html).
