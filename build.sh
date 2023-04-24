#!/bin/bash -f

if [[ -d build ]]
then
    rm -rf build
fi

export SIMPLE_TESTING=FALSE

mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j

