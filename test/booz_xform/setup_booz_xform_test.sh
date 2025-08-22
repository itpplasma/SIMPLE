#!/bin/bash
# Setup script for booz_xform testing
# Downloads test data and prepares comparison environment

set -e  # Exit on error

echo "=== Setting up booz_xform test environment ==="

# Create test directory
TEST_DIR="test/booz_xform"
mkdir -p $TEST_DIR
cd $TEST_DIR

# Download booz_xform file if not present
BOOZ_FILE="boozmn_LandremanPaul2021_QA_lowres.nc"
if [ ! -f "$BOOZ_FILE" ]; then
    echo "Downloading booz_xform test file..."
    wget https://github.com/itpplasma/benchmark_orbit/raw/refs/heads/main/booz_xform/boozmn_LandremanPaul2021_QA_lowres.nc
else
    echo "Booz_xform file already exists: $BOOZ_FILE"
fi

# Download corresponding VMEC file if not present
VMEC_FILE="wout.nc"
if [ ! -f "$VMEC_FILE" ]; then
    echo "Downloading corresponding VMEC file..."
    wget https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc -O wout.nc
else
    echo "VMEC file already exists: $VMEC_FILE"
fi

# Create simple.in for testing
cat > simple.in << EOF
&config
trace_time = 1d-2        ! slowing down time, s
sbeg = 0.5d0             ! starting s (normalized toroidal flux) for particles
send = 0.5d0             ! ending s for uniform sampling
ntestpart = 1            ! number of test particles (1 for coordinate comparison)
netcdffile = 'wout.nc'   ! VMEC file
booz_file = 'boozmn_LandremanPaul2021_QA_lowres.nc'  ! booz_xform file
isw_field_type = 2       ! 2: Internal Boozer conversion (will add 5 for booz_xform)
startmode = 1            ! 1: uniform in s
/
EOF

echo "Setup complete!"
echo "Files created in: $(pwd)"
echo "  - $BOOZ_FILE (booz_xform output)"
echo "  - $VMEC_FILE (VMEC input)"
echo "  - simple.in (test configuration)"