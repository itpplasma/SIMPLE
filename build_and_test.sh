#!/bin/bash -f

if [[ -d build ]]
then
    rm -rf build
fi

export work_dir=$PWD
export FC=gfortran

cd $work_dir

mkdir -p build
cd build
cmake .. -DSIMPLE_TESTING=ON
make -j8

# Execute tests.
echo
echo SIMPLE Tests:
echo
cd ./test/tests
./utility_test.x

cd $work_dir
