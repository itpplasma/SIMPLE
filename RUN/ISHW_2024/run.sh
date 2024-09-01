#!/bin/bash
set -e

cmake --build $CODE/SIMPLE/build -j$(nproc)

pushd vmec
$CODE/SIMPLE/build/simple.x
popd

pushd canflux
$CODE/SIMPLE/build/simple.x
popd

pushd boozer
$CODE/SIMPLE/build/simple.x
popd

pushd meiss
$CODE/SIMPLE/build/simple.x
popd

pushd albert
$CODE/SIMPLE/build/simple.x
popd
