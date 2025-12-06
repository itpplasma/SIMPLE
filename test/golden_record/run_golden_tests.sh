#!/bin/bash

# Run golden record tests for a SIMPLE build
# Usage: run_golden_tests.sh <project_root> <run_dir> <test_data_dir> [single_case]

set -e

PROJECT_ROOT="$1"
RUN_DIR="$2"
TEST_DATA_DIR="$3"
SINGLE_CASE="$4"  # Optional: specific test case to run

if [ -z "$PROJECT_ROOT" ] || [ -z "$RUN_DIR" ] || [ -z "$TEST_DATA_DIR" ]; then
    echo "Usage: $0 <project_root> <run_dir> <test_data_dir> [single_case]"
    exit 1
fi

SIMPLE_EXE="$PROJECT_ROOT/build/simple.x"
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

if [ ! -f "$SIMPLE_EXE" ]; then
    echo "Error: SIMPLE executable not found at $SIMPLE_EXE"
    exit 1
fi

# Create run directory
mkdir -p "$RUN_DIR"

# Find test cases
if [ -n "$SINGLE_CASE" ]; then
    TEST_CASES="$SINGLE_CASE"
    echo "Running single test case: $SINGLE_CASE"
else
    TEST_CASES=$(cd "$SCRIPT_DIR" && find . -name simple.in -exec dirname {} \; | sed 's|^\./||' | grep -v "^\.$" | sort)
    echo "Running all test cases"
fi

echo "Running tests for: $PROJECT_ROOT"
echo "Output directory: $RUN_DIR"
echo "Test cases: $TEST_CASES"
echo "Checking for optional symplectic tokamak regression..."

# Download wout.nc if not present in test data directory
WOUT_FILE="$TEST_DATA_DIR/wout.nc"
if [ ! -f "$WOUT_FILE" ]; then
    echo "Downloading wout.nc..."
    wget -q https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc -O "$WOUT_FILE"
fi

# Download NCSX wout.nc for coils tests
WOUT_NCSX_FILE="$TEST_DATA_DIR/wout_ncsx.nc"
if [ ! -f "$WOUT_NCSX_FILE" ]; then
    echo "Downloading NCSX wout.nc..."
    wget -q https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_c09r00_fixedBoundary_0.5T_vacuum_ns201.nc -O "$WOUT_NCSX_FILE"
fi

# Download and convert NCSX coils file
COILS_STELLOPT="$TEST_DATA_DIR/coils.c09r00"
COILS_SIMPLE="$TEST_DATA_DIR/coils.c09r00.simple"
if [ ! -f "$COILS_SIMPLE" ]; then
    echo "Downloading NCSX coils file..."
    wget -q https://princetonuniversity.github.io/STELLOPT/examples/coils.c09r00 -O "$COILS_STELLOPT"

    echo "Converting coils to simple format..."
    if command -v python3 &> /dev/null; then
        # Try using libneo converter if available
        python3 -m libneo.convert_coils_to_simple "$COILS_STELLOPT" "$COILS_SIMPLE" 2>/dev/null || {
            # Fallback: manual conversion
            echo "Warning: libneo converter not available, using manual conversion"
            python3 - "$COILS_STELLOPT" "$COILS_SIMPLE" <<'EOF'
import sys
with open(sys.argv[1], "r") as f:
    coords = []
    currents = []
    in_filament = False
    for line in f:
        line = line.strip()
        if "begin filament" in line.lower():
            in_filament = True
        elif not in_filament or not line or line.startswith("#"):
            continue
        else:
            parts = line.split()
            if len(parts) >= 4:
                try:
                    x, y, z, I = float(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])
                    coords.append([x, y, z])
                    currents.append(I)
                except ValueError:
                    pass
with open(sys.argv[2], "w") as f:
    f.write(f"{len(coords)}\n")
    for (x, y, z), I in zip(coords, currents):
        f.write(f"{x:.14E}   {y:.14E}   {z:.14E}   {I:.14E}\n")
EOF
        }
    fi
fi

# Run each test case
for test_case in $TEST_CASES; do
    echo "Running test case: $test_case"

    # Create test case directory
    test_dir="$RUN_DIR/$test_case"
    mkdir -p "$test_dir"

    # Copy input file
    cp "$SCRIPT_DIR/$test_case/simple.in" "$test_dir/"

    # Determine which wout and coils files to use based on test case
    if [[ "$test_case" == *"_coils"* ]]; then
        # Coils test cases use NCSX wout and coils files
        ln -sf "$WOUT_NCSX_FILE" "$test_dir/wout.nc"
        ln -sf "$COILS_SIMPLE" "$test_dir/coils.simple"
    elif [[ "$test_case" == *"_ncsx"* ]]; then
        # NCSX test cases use NCSX wout without coils
        ln -sf "$WOUT_NCSX_FILE" "$test_dir/wout.nc"
    else
        # Regular test cases use standard wout
        ln -sf "$WOUT_FILE" "$test_dir/wout.nc"
    fi

    # Run SIMPLE
    cd "$test_dir"
    "$SIMPLE_EXE" > simple.log 2>&1

    if [ $? -ne 0 ]; then
        echo "Error: SIMPLE failed for test case $test_case"
        cat simple.log
        exit 1
    fi

    # Check for expected output files
    if [ ! -f "confined_fraction.dat" ]; then
        echo "Error: confined_fraction.dat not found for test case $test_case"
        exit 1
    fi

    if [ ! -f "times_lost.dat" ]; then
        echo "Error: times_lost.dat not found for test case $test_case"
        exit 1
    fi

    # Run albert diagnostic for albert_coils test cases
    if [[ "$test_case" == "albert_coils" ]]; then
        ALBERT_DIAG_BIN="$PROJECT_ROOT/build/test/tests/field_can/test_field_can_albert_coils_diagnostic.x"
        if [ -x "$ALBERT_DIAG_BIN" ]; then
            echo "Running Albert coordinate diagnostic..."
            "$ALBERT_DIAG_BIN" >> simple.log 2>&1
            if [ $? -ne 0 ]; then
                echo "Warning: Albert diagnostic failed (non-fatal)"
            fi
        fi
        # Also run intermediate values test for detailed regression tracking
        ALBERT_INTER_BIN="$PROJECT_ROOT/build/test/tests/field_can/test_albert_intermediate.x"
        if [ -x "$ALBERT_INTER_BIN" ]; then
            echo "Running Albert intermediate values diagnostic..."
            "$ALBERT_INTER_BIN" >> simple.log 2>&1
            if [ $? -ne 0 ]; then
                echo "Warning: Albert intermediate diagnostic failed (non-fatal)"
            fi
        fi
    fi

    echo "Test case $test_case completed successfully"
done

#
# Optional symplectic tokamak regression (test_sympl_tok.x)
#
SYMP_TOK_BIN="$PROJECT_ROOT/build/test/tests/test_sympl_tok.x"
SYMP_TOK_DIR="$RUN_DIR/sympl_tokamak"
if [ -x "$SYMP_TOK_BIN" ]; then
    echo "Running symplectic tokamak regression via test_sympl_tok.x"
    mkdir -p "$SYMP_TOK_DIR"
    cd "$SYMP_TOK_DIR"
    "$SYMP_TOK_BIN" > sympl_tokamak.log 2>&1
    if [ $? -ne 0 ]; then
        echo "Error: test_sympl_tok.x failed; see $SYMP_TOK_DIR/sympl_tokamak.log"
        exit 1
    fi
    echo "Symplectic tokamak regression completed."
else
    echo "Skipping symplectic tokamak regression (binary not found at $SYMP_TOK_BIN)"
fi

echo "All tests completed successfully"
