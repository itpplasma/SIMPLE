#!/bin/bash
set -e

cmake --build ../build --target simple.x -j$(nproc)

pushd vmec
../../build/simple.x
popd

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
