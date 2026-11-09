#!/bin/bash -f

read -p "Release tag to check out [hash, tag, or branch]: " simple_release_version
echo "Removing old legacy data..."
rm -rf ./simple_*
wget https://raw.githubusercontent.com/hiddenSymmetries/simsopt/912bdcc729d95433518c71df82b8cc45287908b3/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc -O wout.nc

git clone git@github.com:itpplasma/SIMPLE.git simple_$simple_release_version
cd simple_$simple_release_version
git fetch
git checkout ${simple_release_version}
