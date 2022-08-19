#!/bin/bash -f

if [[ -d build ]]
then
    rm -rf build
fi

export SIMPLE_COVERAGE=TRUE

export work_dir=$PWD
export FC=gfortran-8

git submodule update --init --recursive
cd SRC/contrib/pFUnit
#if [[ -d build ]]
#then
#    rm -rf build
#fi
#mkdir build
cd build
#cmake .. -DSKIP_MPI=YES
#make -j16 tests
make -j16 install

cd $work_dir

mkdir -p build
cd build
export PFUNIT_DIR=./SRC/contrib/pFUnit/build/installed/PFUNIT-4.4/cmake
#export PFUNIT_DIR=../../pFUnit/build/installed/PFUNIT-4.4/cmake
cmake .. -DCMAKE_PREFIX_PATH=$PFUNIT_DIR -DCMAKE_INSTALL_PREFIX=$PFUNIT_DIR
make -j16

ctest --output-on-failure
cd SRC/CMakeFiles/SIMPLE.dir/

gcov-9 *.gcno
lcov --gcov-tool gcov-9 --capture --no-recursion --directory . --output-file covered.info
lcov --gcov-tool gcov-9 --capture --no-recursion -i --directory . --output-file uncovered.info
lcov -a covered.info -a uncovered.info --output-file result.info
genhtml --output-directory ../../../COVERAGE result.info

