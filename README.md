# Utilities for and Analysis of feltor simulations

The goal of this collection of notebooks is to provide utilities for
and analysis of the 3d feltor code.

[![LICENSE : MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Installation

The notebooks depend on the `pandas` and `numpy` packages and the `simplesimdb` module. The
latter can be obtained via
```bash
git clone https://github.com/mwiesenberger/simplesimdb
cd simplesimdb
pip install -e .
```

The other dependency is the [Feltor](https://github.com/feltor-dev/feltor) code repository.
Follow the quick-start guide to install.
It is recommended to keep Feltor and this repository next to each other.
If not, you need to set the `FELTOR_PATH` environment variable in order for
the `execute.sh` script to compile and execute the `path/to/feltor/inc/geometries/ds_t` code.

Lastly, we need jupyter to run the notebooks.

## Usage
Run the jupyter notebooks in the repository
```bash
jupyter notebook
```
## Author
Matthias Wiesenberger
