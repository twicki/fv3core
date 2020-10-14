#!/bin/bash

# stop on all erros
set -ex

rundir=$1
if [ ! -d "${rundir}" -o ! -f "${rundir}/diag_table" ] ; then
  echo "Error: must pass an existing rundir as argument"
  exit 1
fi

set +e
grep 'test_case' ${rundir}/input.nml 2>/dev/null 1>/dev/null
if [ $? -eq 0 ] ; then
  is_test_case=1
else
  is_test_case=0
fi

cd "${rundir}"
/bin/rm -rf ./*.out
/bin/rm -rf ./*.log
/bin/rm -rf ./RESTART/*
/bin/rm -rf ./*.nc
/bin/rm -rf ./.gt_cache*
if [ $is_test_case -eq 1 ] ; then
  /bin/rm -rf ./INPUT
  /bin/rm -rf ./grb
  /bin/rm -rf ./co2*.txt
  /bin/rm -rf ./sfc_emiss*.txt
  /bin/rm -rf ./solarconst*.txt
fi
cd -

# done
exit 0
