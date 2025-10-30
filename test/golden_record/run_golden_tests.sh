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

# Download wout.nc if not present in test data directory
WOUT_FILE="$TEST_DATA_DIR/wout.nc"
if [ ! -f "$WOUT_FILE" ]; then
    echo "Downloading wout.nc..."
    wget -q https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc -O "$WOUT_FILE"
fi

# Run each test case
for test_case in $TEST_CASES; do
    echo "Running test case: $test_case"
    
    # Create test case directory
    test_dir="$RUN_DIR/$test_case"
    mkdir -p "$test_dir"
    
    # Copy input file
    cp "$SCRIPT_DIR/$test_case/simple.in" "$test_dir/"
    
    # Create symlink to wout.nc
    ln -sf "$WOUT_FILE" "$test_dir/wout.nc"
    
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
    
    echo "Test case $test_case completed successfully"
done

echo "All tests completed successfully"