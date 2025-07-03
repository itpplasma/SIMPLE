#!/usr/bin/env python3
"""
Convert VMEC wout file to GVEC state file for testing.
This script uses GVEC's Python interface to read VMEC and save as GVEC format.
"""

import sys
import os

try:
    # Add GVEC Python path if needed
    gvec_path = os.path.join(os.path.dirname(__file__), 
                            "../../../build/_deps/gvec-build/python")
    if os.path.exists(gvec_path):
        sys.path.insert(0, gvec_path)
    
    import gvec
    
    # Read VMEC file
    print("Reading VMEC file: wout.nc")
    gvec.read_vmec("wout.nc", file_format=0)  # 0 for NetCDF
    
    # Initialize GVEC from VMEC data
    print("Initializing GVEC from VMEC data")
    gvec.init_vmec()
    
    # Save state
    output_file = "gvec_from_vmec.dat"
    print(f"Saving GVEC state to: {output_file}")
    gvec.save(output_file)
    
    print("Conversion successful!")
    sys.exit(0)
    
except ImportError as e:
    print(f"ERROR: Could not import gvec module: {e}")
    print("GVEC Python module not available")
    sys.exit(1)
except Exception as e:
    print(f"ERROR: {e}")
    sys.exit(1)