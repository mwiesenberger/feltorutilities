#!/bin/bash

: ${FELTOR_PATH:="../../feltor"}
# If feltor is not here then change the FELTOR_PATH enviromnent variable
# export FELTOR_PATH="path/to/feltor"

make -C $FELTOR_PATH/src/feltor/ feltor

rm -f $2
$FELTOR_PATH/src/feltor/feltor $1 $2
