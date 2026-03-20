#!/bin/bash
set -euo pipefail

SIMPLE_X="$(cd "$(dirname "$0")/../.." && pwd)/build/simple.x"
EXAMPLE_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="/tmp/simple_classification"

if [ ! -x "$SIMPLE_X" ]; then
    echo "ERROR: simple.x not found at $SIMPLE_X"
    echo "Build with: make -C $(dirname "$SIMPLE_X")/.."
    exit 1
fi

# Equilibrium files shipped in data/ directory
DATA_DIR="$EXAMPLE_DIR/data"
QI_WOUT="$DATA_DIR/wout_23_1900_fix_bdry.nc"
QH_WOUT="$DATA_DIR/wout_qh_8_7.nc"
QA_WOUT="$DATA_DIR/wout_henneberg_qa.nc"

for f in "$QI_WOUT" "$QH_WOUT" "$QA_WOUT"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Equilibrium file not found: $f"
        echo "These should be in examples/classification/data/"
        exit 1
    fi
done

mkdir -p "$BASE_DIR"

# Copy equilibrium files to base directory (avoids symlink issues)
cp -u "$QI_WOUT" "$BASE_DIR/wout_23_1900_fix_bdry.nc"
cp -u "$QH_WOUT" "$BASE_DIR/wout_qh_8_7.nc"
cp -u "$QA_WOUT" "$BASE_DIR/wout_henneberg_qa.nc"

run_case() {
    local case_name="$1"
    local dir="$BASE_DIR/$case_name"

    mkdir -p "$dir"
    cp "$EXAMPLE_DIR/configs/${case_name}.in" "$dir/simple.in"

    # Link all wout files into the run directory
    for wout in "$BASE_DIR"/wout_*.nc; do
        ln -sf "$wout" "$dir/"
    done

    echo "=== Running $case_name ==="
    local t_start
    t_start=$(date +%s)
    (cd "$dir" && "$SIMPLE_X")
    local t_end
    t_end=$(date +%s)
    local elapsed=$((t_end - t_start))
    echo "=== $case_name done in ${elapsed}s ==="
}

# Determine which cases to run
if [ $# -gt 0 ]; then
    # Run only specified cases
    for case_name in "$@"; do
        run_case "$case_name"
    done
else
    # Run all cases: loss runs first (Fig 5-7), then Fig 8, then volume
    LOSS_CASES="qi_s03 qi_s06 qh_s03 qh_s06 qa_s03 qa_s06"
    FIG8_CASES="qi_fig8 qh_fig8 qa_fig8"
    VOLUME_CASES="qh_volume"

    for case_name in $LOSS_CASES $FIG8_CASES $VOLUME_CASES; do
        run_case "$case_name"
    done
fi

echo ""
echo "All runs complete. Results in $BASE_DIR"
