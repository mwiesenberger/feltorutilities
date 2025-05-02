# Utilities for feltor simulations

The goal of this package is to provide utilities for the setup and run
 of 3d feltor simulations.

It consists of
 - a Python module `feltorutilities` containing parameter functions and utilities to
convert from dimensional to dimensionless parameters and back.
 - Jupyter notebooks guiding the generation of suitable input parameters for both the feltor and the thermal feltor codes.


[![LICENSE : MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Installation
### Python module
The python package `feltorutilities` can be locally installed via
```bash
git clone https://github.com/mwiesenberger/feltorutilities
cd feltorutilities
pip install -e .
```

### Jupyter notebooks
The notebooks depend on the `pandas` and `numpy` packages and the `simplesimdb` and `pyfeltor` module.
The other dependency is the [Feltor](https://github.com/feltor-dev/feltor) code repository.
Follow the quick-start guide to install.
It is recommended to keep Feltor and this repository next to each other.
If not, you will need to set the `FELTOR_PATH` environment variable in order for
`feltor.sh`, `job.sh` and `feltordiag.sh` to work.

Lastly, we need jupyter to run the notebooks.

## Usage
### A plasma formulary

First, a list of plasma formulas is provided. These functions can be called as
```python
import feltorutilities as fp
fp.rho_s( B_0, m_i, T_e)
# or
parameters = {"B_0" : 0.7, "m_i" : fp.deuteron_mass, "T_e" : 0.07}
fp.rho_s( **parameters)
```
To get a list of available functions type

```python
function_names = fp.quantities()
```

A utility function is provided that lets you compute a collection of quantities
from given parameters at once, e.g.

```python
parameters = {"B_0" : 0.7, "m_i" : fp.deuteron_mass, "T_e" : 0.07}
values = parameters2quantities( parameters, ["rho_s", "c_s", "omega_0_inv"])
```

Second we provide a utility function that inverts given numerical (dimensionless) parameters
back to the corresponding physical (dimensional) ones. This is useful in analysing
simulation data where only the numerical values are stored.
For example:
```python
physical = {"R" : 0.6}
numerical = {"mu" : -2.72e-4,"tau":1,"beta":4.11e-4,"resistivity" : 3.81e-5,"R_0" : 91.94}
numerical2physical( numerical, physical)
# now physical contains R, m_i, T_e, n_0, B_0, T_i)
```

### Feltor simulation setup
This package contains two folders `feltor` and `thermal` that are intended to serve as templates to start your own project for the corresponding feltor codes.
In order to setup a Feltor simulation you need to
- First, decide if you want to run a isothermal, 2-species Feltor simulation  or a thermal, multispecies simulation. In the first case copy the `feltor` folder to a new repository, in the second one use the `thermal` folder
- If necessary adapt the `FELTOR_PATH` variable to the feltor C++ repository as described above in all `*.sh` files. The explanation for the structure of these files can be found in the `simplesimdb` package.
- Run the jupyter notebook in the folder and follow the instructions in the file
- Once you generated satisfactory input parameters, either manually copy them or load them in the `generate_data.py` script. The idea for this script is to use the `simplesimdb` package to automatically generate submit scripts on a HPC cluster such that the simulation is automatically restarted if necessary and/or parameter studies are possible.

## Author
Matthias Wiesenberger
