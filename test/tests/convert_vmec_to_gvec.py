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
    
    # First, let's see what's available in the GVEC module
    print("GVEC module attributes:")
    attrs = [attr for attr in dir(gvec) if not attr.startswith('_')]
    for attr in attrs:
        print(f"  {attr}")
    
    # Check if there's a specific VMEC-related module or function
    if hasattr(gvec, 'vmec'):
        print("\nGVEC.vmec attributes:")
        vmec_attrs = [attr for attr in dir(gvec.vmec) if not attr.startswith('_')]
        for attr in vmec_attrs:
            print(f"  {attr}")
    
    # Check Run class
    if hasattr(gvec, 'Run'):
        print("\nGVEC.Run attributes:")
        run_attrs = [attr for attr in dir(gvec.Run) if not attr.startswith('_')]
        for attr in run_attrs:
            print(f"  {attr}")
    
    # Check State class
    if hasattr(gvec, 'State'):
        print("\nGVEC.State attributes:")
        state_attrs = [attr for attr in dir(gvec.State) if not attr.startswith('_')]
        for attr in state_attrs:
            print(f"  {attr}")
    
    # Check if there's a run function that might help
    if hasattr(gvec, 'run'):
        print("\nChecking gvec.run function...")
        print(f"  Type: {type(gvec.run)}")
        if hasattr(gvec.run, '__doc__'):
            print(f"  Docstring: {gvec.run.__doc__}")
    
    # Try to create a GVEC run from VMEC
    vmec_file = "wout.nc"
    output_file = "wout_gvec.dat"
    
    print(f"\nAttempting to convert VMEC file: {vmec_file}")
    
    # Try using the Run class
    if hasattr(gvec, 'Run'):
        try:
            print("Trying to create GVEC Run from VMEC...")
            # Try different initialization methods
            run = gvec.Run()  # Create empty run
            if hasattr(run, 'from_vmec'):
                run.from_vmec(vmec_file)
            elif hasattr(run, 'load_vmec'):
                run.load_vmec(vmec_file)
            elif hasattr(run, 'read_vmec'):
                run.read_vmec(vmec_file)
            else:
                print("  Run object methods:")
                for attr in dir(run):
                    if not attr.startswith('_'):
                        print(f"    {attr}")
        except Exception as e:
            print(f"  Failed with Run class: {e}")
    
    # Try using load_state if available
    if hasattr(gvec, 'load_state'):
        print("\nChecking load_state function...")
        print(f"  Type: {type(gvec.load_state)}")
        if hasattr(gvec.load_state, '__doc__'):
            print(f"  Docstring: {gvec.load_state.__doc__}")
    
    print("\nERROR: Could not determine correct GVEC API for VMEC conversion")
    print("Please check GVEC documentation for the correct Python API")
    sys.exit(1)
    
    print("Conversion successful!")
    sys.exit(0)
    
except ImportError as e:
    print(f"ERROR: Could not import gvec module: {e}")
    print("GVEC Python module not available")
    sys.exit(1)
except Exception as e:
    print(f"ERROR: {e}")
    sys.exit(1)