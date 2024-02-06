#!/bin/bash -f

if [[ -d build ]]
then
    rm -rf build
fi

mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release .. -DSIMPLE_TESTING=OFF
make -j8
