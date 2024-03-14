#!/bin/bash -f

read -p "Release tag to check out: " simple_release_version

git clone git@github.com:itpplasma/SIMPLE.git simple_$simple_release_version
cd simple_$simple_release_version
git fetch --all --tags --prune
git checkout tags/${simple_release_version}