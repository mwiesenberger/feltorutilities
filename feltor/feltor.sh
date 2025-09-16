#!/bin/bash

# We suggest a serial feltor version to run with the setup
: ${FELTOR_PATH:="../../feltor/build/cpu"}
# If feltor is not here then change the FELTOR_PATH enviromnent variable
# export FELTOR_PATH="path/to/feltor/build/X"
#
# In order to compile use
# cd path/to/feltor
# cmake --preset cpu
# cmake --build build/cpu --target feltor_feltor
rm -f $2
$FELTOR_PATH/src/feltor/feltor $1 $2
