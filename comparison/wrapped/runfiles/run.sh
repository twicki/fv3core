#!/bin/bash

set -ex

ulimit -s
mpirun -l -np 6 python model.py 2>stderr.log | tee stdout.log

grep Termination 2>/dev/null 1>/dev/null stdout.log
if [ $? -ne 0 ] ; then
    echo "Error: problem with this model run"
    exit 1
fi

exit 0
