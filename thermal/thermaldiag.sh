#!/bin/bash

# We suggest a serial feltor version to run with diagnostics
# (runtime depends on harddrive speed alone)
: ${FELTOR_PATH:="../../feltor/build/cpu"}
# If feltor is not here then change the FELTOR_PATH enviromnent variable
# export FELTOR_PATH="path/to/feltor/build/X"
#
# In order to compile use
# cd path/to/feltor
# cmake --preset cpu
# cmake --build build/cpu --target thermal_thermaldiag

input=$(echo $2 | sed -e 's/diag/thermal3d/')

$FELTOR_PATH/src/thermal/thermaldiag "config.json" $input $2

