# Calkulate

Calkulate is a Python package for finding total alkalinity from titration data using [PyCO2SYS](https://pyco2sys.rtfd.io).

## Installation

    pip install calkulate

## Usage

If the data for each individual titration is in its own text file and you have a spreadsheet containing the metadata for each titration — [all formatted as expected]() — then all you need to do with Calkulate is:

```python
import calkulate as calk
data = calk.read_csv("path/to/metadata_file.csv").calkulate()
```

`data` is then a pandas DataFrame based on the metadata file you provided but with some extra columns added such as `data.alkalinity`, which contains the fully calibrated total alkalinity for each sample.

Other read-in functions are also available (e.g. for [Excel spreadsheets]() and [VINDTA .dbs]() files).

## About

Calkulate is being developed by [Dr Matthew Humphreys](https://mvdh.xyz) at the Royal Netherlands Institute for Sea Research ([NIOZ, Texel, the Netherlands](https://www.nioz.nl/en)).

## Citation

If you use Calkulate in your work, please cite it as:

!!! info "Calkulate citation"
    Humphreys, M. P. and Matthews, R. S. (2020).  Calkulate: total alkalinity from titration data in Python.  *Zenodo.*  [doi:10.5281/zenodo.2634304](https://doi.org/10.5281/zenodo.2634304).

Please specify which version you are using.  To find this:

```python
import calkulate as calk
calk.say_hello()
```

## License

Calkulate is licensed under the [GNU General Public License version 3 (GPLv3)](https://www.gnu.org/licenses/gpl-3.0.en.html).
