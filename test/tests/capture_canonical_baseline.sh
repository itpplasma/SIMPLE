#!/bin/bash
# Script to capture canonical coordinates baseline for regression testing

echo "Capturing canonical coordinates baseline..."

# Ensure we have the VMEC file
if [ ! -f "wout.nc" ]; then
    echo "Downloading VMEC test file..."
    wget -q https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc -O wout.nc
fi

# Ensure we have the input file configured for canonical coordinates
cat > simple.in << EOF
&config
trace_time = 1d-6        ! slowing down time, s (very short for debug)
sbeg = 0.3d0             ! starting s (normalzed toroidal flux) for particles.
ntestpart = 2            ! number of test particles (just 2 for debug)
netcdffile = 'wout.nc'   ! name of VMEC file in NETCDF format
isw_field_type = 0       ! -1: Testing, 0: Canonical, 1: VMEC, 2: Boozer, 3: Meiss, 4: Albert
/
EOF

echo "Running simple.x to capture canonical coordinates debug output..."
./build/simple.x | grep "DEBUG:" > canonical_baseline.txt 2>&1

echo "Baseline captured in canonical_baseline.txt:"
cat canonical_baseline.txt

echo "Done. Use this baseline for regression testing."