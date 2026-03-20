#!/bin/bash
# Run all classification cases from the JPP 2020 paper.
# Equilibrium files are downloaded automatically on first run.
set -euo pipefail

SIMPLE_X="$(cd "$(dirname "$0")/../.." && pwd)/build/simple.x"
EXAMPLE_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${1:-/tmp/simple_classification}"

if [ ! -x "$SIMPLE_X" ]; then
    echo "ERROR: simple.x not found. Build with: make"
    exit 1
fi

mkdir -p "$BASE_DIR"

# Download equilibrium files if not already present
download_wout() {
    local name="$1" url="$2" dest="$BASE_DIR/$name"
    if [ ! -f "$dest" ]; then
        echo "Downloading $name ..."
        wget -q -O "$dest" "$url"
    fi
}

download_wout wout_qh_8_7.nc \
    "https://github.com/itpplasma/proxima-simple-classification/raw/main/data/wout_qh_8_7.nc"
download_wout wout_23_1900_fix_bdry.nc \
    "https://github.com/itpplasma/proxima-simple-classification/raw/main/data/wout_23_1900_fix_bdry.nc"
download_wout wout_henneberg_qa.nc \
    "https://github.com/itpplasma/proxima-simple-classification/raw/main/data/wout_henneberg_qa.nc"

run_case() {
    local case_name="$1"
    local dir="$BASE_DIR/$case_name"

    mkdir -p "$dir"
    cp "$EXAMPLE_DIR/configs/${case_name}.in" "$dir/simple.in"

    for wout in "$BASE_DIR"/wout_*.nc; do
        ln -sf "$wout" "$dir/"
    done

    echo "=== Running $case_name ==="
    local t0; t0=$(date +%s)
    (cd "$dir" && "$SIMPLE_X")
    echo "=== $case_name done in $(($(date +%s) - t0))s ==="
}

if [ $# -gt 1 ]; then
    shift  # skip base_dir
    for case_name in "$@"; do run_case "$case_name"; done
else
    for case_name in \
        qi_s03 qi_s06 qh_s03 qh_s06 qa_s03 qa_s06 \
        qi_fig8 qh_fig8 qa_fig8 \
        qi_volume qh_volume qa_volume; do
        run_case "$case_name"
    done
fi

echo "All runs complete. Results in $BASE_DIR"
