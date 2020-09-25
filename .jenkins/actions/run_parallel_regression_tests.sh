#!/bin/bash
set -e -x
BACKEND=$1
EXPNAME=$2
ARGS="-vv -s -rsx --backend=${BACKEND} --which_modules=CubedToLatLon,DynCore,FVSubgridZ,HaloUpdate,HaloUpdate_2,HaloVectorUpdate,MPPBoundaryAdjust,MPPUpdateDomains,Tracer2D1L --junitxml=/.jenkins/parallel_test_results.xml"
export EXPERIMENT=${EXPNAME}

# Set the host data location
export TEST_DATA_HOST="${TEST_DATA_DIR}/${EXPNAME}/"

# sync the test data
make get_test_data

# The default of this set to 1 causes a segfault
make run_tests_parallel TEST_ARGS="${ARGS}"

ARGS_FVDYN="-vv -s -rsx --backend=${BACKEND} --which_modules=FVDynamics --junitxml=/.jenkins/parallel_test_results.xml"
make run_tests_parallel TEST_ARGS="${ARGS_FVDYN}"

echo `ls -lh ${TEST_DATA_HOST}/*.txt`
echo `cat ${TEST_DATA_HOST}/regression*.txt`
