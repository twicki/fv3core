#!/bin/bash
set -e -x
echo $LD_LIBRARY_PATH
echo $PATH
echo `ls -lh /project/s1053/install/modulefiles/`
echo `ls -lh /project/s1053/install/modulefiles/gcloud`
module load daint-gpu
module add "/project/s1053/install/modulefiles/"
module load gcloud
echo `which gsutil`
