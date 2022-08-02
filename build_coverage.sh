#!/bin/bash -f

if [[ -d build ]]
then
    rm -rf build
fi

export SIMPLE_COVERAGE=TRUE

mkdir -p build
cd build
cmake .. -DCMAKE_PREFIX_PATH=$PFUNIT_DIR
make -j

ctest --output-on-failure
cd SRC/CMakeFiles/SIMPLE.dir/

gcov-9 *.gcno
lcov --gcov-tool gcov-9 --capture --no-recursion --directory . --output-file covered.info
lcov --gcov-tool gcov-9 --capture --no-recursion -i --directory . --output-file uncovered.info
lcov -a covered.info -a uncovered.info --output-file result.info
genhtml --output-directory ../../../COVERAGE result.info

