# Calkulate

Calkulate is a Python package for finding total alkalinity from titration data.

## Installation

This documentation site is for version 3, which is in development.  For now, install with: 

    pip install git+https://github.com/mvdh7/calkulate.git@v3

For [version 2](https://calkulate.readthedocs.io/en/latest/), and for version 3 once it's ready, install with:

    pip install calkulate

## Get started

  1.  Read how to import and work with data from [a single titration](../titrations).
  2.  See the additional tools for efficiently investigating [datasets of multiple titrations](../datasets).

## About

Calkulate is being developed by [Dr Matthew Humphreys](https://mvdh.xyz) at the Royal Netherlands Institute for Sea Research ([NIOZ, Texel, the Netherlands](https://www.nioz.nl/en)).

## Citation

If you use Calkulate in your work, please cite it as:

!!! info "Calkulate citation"
    Humphreys, M. P. and Matthews, R. S. (2020).  Calkulate: total alkalinity from titration data in Python.  *Zenodo.*  [doi:10.5281/zenodo.2634304](https://doi.org/10.5281/zenodo.2634304).

Please specify which version you are using.  To find this:

    :::python
    import calkulate as calk
    calk.say_hello()

## License

Calkulate is licensed under the [GNU General Public License version 3 (GPLv3)](https://www.gnu.org/licenses/gpl-3.0.en.html).
