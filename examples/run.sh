#!/bin/bash

cmake --build ../build --target simple.x -j$(nproc)

pushd canflux
../../build/simple.x
popd

pushd boozer
../../build/simple.x
popd

pushd meiss
../../build/simple.x
popd

pushd albert
../../build/simple.x
popd
