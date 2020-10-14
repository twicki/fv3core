#!/bin/bash

set -ex

# make the environment sane
export DEBIAN_FRONTEND=noninteractive
apt-get update
apt-get -y install bc vim netcdf-bin cdo nco tkdiff
ln -fs /usr/share/zoneinfo/America/Los_Angeles /etc/localtime
dpkg-reconfigure --frontend noninteractive tzdata

# required to install packages in develop mode
python setup.py develop

# we need newer version of Serialbox
git clone -b v2.6.1 --depth 1 https://github.com/GridTools/serialbox.git /tmp/serialbox
cd /tmp/serialbox
cmake -B build -S /tmp/serialbox -DSERIALBOX_USE_NETCDF=ON -DSERIALBOX_TESTING=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/serialbox
cmake --build build/ -j $(nproc) --target install
cd -
rm -rf /tmp/serialbox

# all done
exit 0

