#!/bin/bash
# Run booz_xform comparison test in build directory
# This script is called by CTest

set -e  # Exit on error

echo "=== Running booz_xform comparison test ==="
echo "Working directory: $(pwd)"

# Download test files if not present
BOOZ_FILE="boozmn_LandremanPaul2021_QA_lowres.nc"
VMEC_FILE="wout.nc"

if [ ! -f "$BOOZ_FILE" ]; then
    echo "Downloading booz_xform test file..."
    wget -q https://github.com/itpplasma/benchmark_orbit/raw/refs/heads/main/booz_xform/boozmn_LandremanPaul2021_QA_lowres.nc
fi

if [ ! -f "$VMEC_FILE" ]; then
    echo "Downloading VMEC test file..."
    wget -q https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc -O wout.nc
fi

# Create simple.in for testing
cat > simple.in << EOF
&config
trace_time = 1d-2
sbeg(1) = 0.5d0
ntestpart = 1
netcdffile = 'wout.nc'
isw_field_type = 2       ! Internal Boozer
/
EOF

# Get the source directory for the Python script
if [ -z "$CMAKE_CURRENT_SOURCE_DIR" ]; then
    # If not set, assume we're in build/test/tests
    SCRIPT_DIR="../../../test/booz_xform"
else
    SCRIPT_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../booz_xform"
fi

# Run the Python comparison script
echo "Running comparison..."
python3 ${SCRIPT_DIR}/compare_boozer_coords.py

# Check if output was created
if [ -f "boozer_comparison.png" ]; then
    echo "Test completed successfully. Output: boozer_comparison.png"
    exit 0
else
    echo "Error: Expected output file not created"
    exit 1
fi