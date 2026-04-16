#!/bin/bash
# Run a quick SIMPLE example: download test data, trace orbits, print result.
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
SIMPLE_X="${1:-$REPO_ROOT/build/simple.x}"
WOUT_URL="https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"

if [ ! -x "$SIMPLE_X" ]; then
    echo "Executable not found: $SIMPLE_X"
    echo "Build first with 'make' from the repository root, or pass the path as argument."
    exit 1
fi

RUN_DIR="/tmp/simple_example"
mkdir -p "$RUN_DIR"
cp "$SCRIPT_DIR/simple.in" "$RUN_DIR/"

if [ ! -f "$RUN_DIR/wout.nc" ]; then
    echo "Downloading test VMEC file..."
    wget -q "$WOUT_URL" -O "$RUN_DIR/wout.nc"
fi

echo "Running SIMPLE in $RUN_DIR ..."
(cd "$RUN_DIR" && "$SIMPLE_X")

echo -e "\nConfined fraction over time:"
column -t "$RUN_DIR/confined_fraction.dat"
