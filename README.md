# Utilities for and Analysis of feltor simulations

The goal of this collection of notebooks is to provide utilities for the setup
and analysis of the 3d feltor code.

The python module `feltorutilities` contains parameter functions and methods to
convert from physical to numerical parameters and back.

The `generate_data.py` script is used to generate the simulation data for a
resistivity study on COMPASS.  The input parameters were setup in the notebook
`COMPASS-setup.ipynb`


[![LICENSE : MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Installation

The notebooks depend on the `pandas` and `numpy` packages and the `simplesimdb` and `pyfeltor` module. The
package can be locally installed via
```bash
git clone https://github.com/mwiesenberger/feltorutilities
cd feltorutilities
pip install -e .
```

The other dependency is the [Feltor](https://github.com/feltor-dev/feltor) code repository.
Follow the quick-start guide to install.
It is recommended to keep Feltor and this repository next to each other.
If not, you need to set the `FELTOR_PATH` environment variable in order for
`feltor.sh` and `feltordiag.sh` to work.

Lastly, we need jupyter to run the notebooks.

## Usage
Run the jupyter notebooks in the repository
```bash
jupyter-lab
```
## Author
Matthias Wiesenberger
