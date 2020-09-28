#!/bin/bash
set -e -x
BACKEND=$1
EXPNAME=$2
echo 'TESTING SUBSET'
ARGS="-vv -s -rsx --backend=${BACKEND} --print_failures --which_modules=CubedToLatLon,DynCore,FVSubgridZ,HaloUpdate,HaloUpdate-2,HaloVectorUpdate,MPPBoundaryAdjust,MPPUpdateDomains,Tracer2D1L --junitxml=/.jenkins/parallel_test_results.xml"
export EXPERIMENT=${EXPNAME}

# Set the host data location
export TEST_DATA_HOST="${TEST_DATA_DIR}/${EXPNAME}/"

# sync the test data 
make get_test_data

# The default of this set to 1 causes a segfault
#make run_tests_parallel TEST_ARGS="${ARGS}"
#echo 'TESTING FVDynamics'
#ARGS_FVDYN="-vv -s -rsx --backend=${BACKEND}  --print_failures  --which_modules=FVDynamics --junitxml=/.jenkins/parallel_test_results.xml"
#make run_tests_parallel TEST_ARGS="${ARGS_FVDYN}"
if [ ${EXPNAME} == "c12_6ranks_standard" ]; then
    make run_tests_parallel TEST_ARGS="${ARGS_FVDYN} -vv -s -rsx --backend=${BACKEND}  --print_failures  --which_modules=FVDynamics --junitxml=/.jenkins/parallel_test_results.xml "
else
for COUNT in 1 2 3 4 5 6
do
    echo "THIS TIME", ${COUNT}
    make run_tests_parallel TEST_ARGS=" -vv -s -rsx --backend=${BACKEND} --which_modules=CubedToLatLon,DynCore,FVSubgridZ,HaloUpdate,HaloUpdate-2,HaloVectorUpdate,MPPBoundaryAdjust,MPPUpdateDomains,Tracer2D1L  --print_failures --faulure_stride=79  --junitxml=/.jenkins/parallel_test_results.xml "
    make run_tests_parallel TEST_ARGS=" -vv -s -rsx --backend=${BACKEND} --which_modules=FVDynamics  --print_failures --faulure_stride=79  --junitxml=/.jenkins/parallel_test_results.xml "
done
fi
echo `ls -lh ${TEST_DATA_HOST}/*.txt`
echo `cat ${TEST_DATA_HOST}/regression*.txt`
