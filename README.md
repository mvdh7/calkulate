# Calkulate

![Tests](https://github.com/mvdh7/calkulate/workflows/Tests/badge.svg)
[![PyPI version](https://img.shields.io/pypi/v/calkulate.svg?style=popout)](https://pypi.org/project/calkulate/)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/calkulate.svg)](https://anaconda.org/conda-forge/calkulate)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.2634304-informational)](https://doi.org/10.5281/zenodo.2634304)
[![Docs](https://readthedocs.org/projects/calkulate/badge/?version=latest&style=flat)](https://calkulate.readthedocs.io/en/latest/)
[![Coverage](https://github.com/mvdh7/calkulate/blob/main/.misc/coverage.svg)](https://github.com/mvdh7/calkulate/blob/main/.misc/coverage.txt)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

Calkulate is a Python package for finding total alkalinity from titration data using [PyCO2SYS](https://PyCO2SYS.rtfd.io).

## Installation

    pip install calkulate
    conda install conda-forge :: calkulate

## Use

If the data for each individual titration is in its own text file and you have a spreadsheet containing the metadata for each titration on separate rows — all formatted as expected — then all you need to do with Calkulate is:

```python
import calkulate as calk
data = calk.read_csv("path/to/metadata_file.csv").calkulate()
data.alkalinity  # <== here are your alkalinity results
```

For more detail, see [the online documentation](https://mvdh.xyz/calkulate/).

## About

Calkulate is being developed primarily by [Dr Matthew P. Humphreys](https://www.nioz.nl/en/about/organisation/staff/matthew-humphreys) at the Royal Netherlands Institute for Sea Research ([NIOZ, Texel](https://www.nioz.nl/en)).

## Citation

If you use Calkulate in your work, please cite it as:

> Humphreys, M. P. and Matthews, R. S. (2025).  Calkulate: total alkalinity from titration data in Python.  *Zenodo.*  [doi:10.5281/zenodo.2634304](https://doi.org/10.5281/zenodo.2634304).

Please report which version you are using.  To find this out:

```python
import calkulate as calk
calk.hello()
```

## License

Calkulate is licensed under the [GNU General Public License version 3 (GPLv3)](https://www.gnu.org/licenses/gpl-3.0.en.html).
