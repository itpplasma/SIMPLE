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
sbeg = 0.5d0
send = 0.5d0
ntestpart = 1
netcdffile = 'wout.nc'
isw_field_type = 2       ! Internal Boozer
startmode = 1
/
EOF

# Run the Python comparison script
echo "Running comparison..."
python3 ${CMAKE_CURRENT_SOURCE_DIR}/../booz_xform/compare_boozer_coords.py

# Check if output was created
if [ -f "boozer_comparison.png" ]; then
    echo "Test completed successfully. Output: boozer_comparison.png"
    exit 0
else
    echo "Error: Expected output file not created"
    exit 1
fi