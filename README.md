# Calkulate

![Tests](https://github.com/mvdh7/calkulate/workflows/Tests/badge.svg)
[![PyPI version](https://img.shields.io/pypi/v/calkulate.svg?style=popout)](https://pypi.org/project/calkulate/)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.2634304-informational)](https://doi.org/10.5281/zenodo.2634304)
[![Docs](https://readthedocs.org/projects/calkulate/badge/?version=latest&style=flat)](https://calkulate.readthedocs.io/en/latest/)
[![Coverage](https://github.com/mvdh7/calkulate/blob/main/.misc/coverage.svg)](https://github.com/mvdh7/calkulate/blob/main/.misc/coverage.txt)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

Calkulate is a Python package for finding total alkalinity from titration data using [PyCO2SYS](https://PyCO2SYS.rtfd.io).

## Installation

    pip install calkulate

## Use

If the data for each individual titration is in its own text file and you have a spreadsheet containing the metadata for each titration on separate rows — all formatted as expected — then all you need to do with Calkulate is:

```python
import calkulate as calk
data = calk.read_csv("path/to/metadata_file.csv").calkulate()
data.alkalinity  # <== here are your alkalinity results
```

For more detail, see [the online documentation](https://calkulate.readthedocs.io/en/latest/).

## About

Calkulate is being developed primarily by [Dr Matthew P. Humphreys](https://mvdh.xyz) at the Royal Netherlands Institute for Sea Research ([NIOZ, Texel, the Netherlands](https://www.nioz.nl/en)).

While its results should be reliable, the package is still a work-in-progress intended primarily for the developers' own use, and not everything is documented.  You are therefore strongly encouraged to get in touch if you're using it.  Contributions are welcome!

## Citation

If you use Calkulate in your work, please cite it as:

> Humphreys, M. P. and Matthews, R. S. (2022).  Calkulate: total alkalinity from titration data in Python.  *Zenodo.*  [doi:10.5281/zenodo.2634304](https://doi.org/10.5281/zenodo.2634304).

Please report which version you are using.  To find this out:

```python
import calkulate as calk
calk.hello()
```

## License

Calkulate is licensed under the [GNU General Public License version 3 (GPLv3)](https://www.gnu.org/licenses/gpl-3.0.en.html).
