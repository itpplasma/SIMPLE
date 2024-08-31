#!/bin/bash
set -e

cmake --build ../build -j$(nproc)

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
