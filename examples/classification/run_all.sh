#!/bin/bash
set -euo pipefail

SIMPLE_X="$(cd "$(dirname "$0")/../.." && pwd)/build/simple.x"
EXAMPLE_DIR="$(cd "$(dirname "$0")" && pwd)"
WOUT_URL="https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc"

BASE_DIR="/tmp/simple_classification"
mkdir -p "$BASE_DIR"

# Download wout if needed
WOUT="$BASE_DIR/wout_LandremanPaul2021_QH_reactorScale_lowres_reference.nc"
if [ ! -f "$WOUT" ]; then
    echo "Downloading VMEC equilibrium file..."
    wget -q -O "$WOUT" "$WOUT_URL"
fi

for case in s03 s06 volume; do
    DIR="$BASE_DIR/$case"
    mkdir -p "$DIR"
    cp "$EXAMPLE_DIR/simple_${case}.in" "$DIR/simple.in"
    ln -sf "$WOUT" "$DIR/"
    echo "=== Running $case ==="
    cd "$DIR"
    "$SIMPLE_X"
    echo "=== $case done ==="
done

echo "All runs complete. Results in $BASE_DIR"
