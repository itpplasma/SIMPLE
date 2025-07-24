#!/usr/bin/env python3
"""
Compare magnetic field evaluation between SIMPLE's internal Boozer conversion 
and BOOZXFORM external file at the same coordinates.

This script:
1. Runs SIMPLE with isw_field_type=2 (internal Boozer) and saves field data
2. Runs SIMPLE with isw_field_type=5 (BOOZXFORM) and saves field data  
3. Compares the field evaluations at the same coordinates
"""
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import sys
import os
import subprocess
import tempfile

# Add parent directory to path for pysimple
script_dir = os.path.dirname(os.path.abspath(__file__))
if 'build' in script_dir:
    # We're running from build directory
    sys.path.insert(0, os.getcwd())  # Current dir should have pysimple
else:
    # We're running from source directory
    build_dir = os.path.join(script_dir, '../../build')
    if os.path.exists(build_dir):
        sys.path.insert(0, build_dir)

try:
    import pysimple
except ImportError as e:
    print(f"Error: pysimple not found. Please build SIMPLE first.")
    print(f"Python path: {sys.path}")
    print(f"Import error: {e}")
    sys.exit(1)


def create_test_config(field_type, boozxform_file=None):
    """Create simple.in configuration for given field type"""
    config = f"""&config
    isw_field_type = {field_type}
    netcdffile = 'wout.nc'
    ntestpart = 1
    nper = 1
    trace_time = 1d-4
    ntimstep = 10
    generate_start_only = .true.
    startmode = 5
    num_surf = 1
    sbeg(1) = 0.5
"""
    if boozxform_file:
        config += f"    boozxform_file = '{boozxform_file}'\n"
    
    config += "/\n"
    return config


def evaluate_field_at_coordinates(field_type, coords, boozxform_file=None):
    """Evaluate magnetic field at given coordinates using specified field type"""
    
    # Create temporary configuration
    config_content = create_test_config(field_type, boozxform_file)
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.in', delete=False) as f:
        f.write(config_content)
        config_file = f.name
    
    try:
        # Read configuration
        pysimple.params.read_config(config_file)
        
        # Initialize tracer
        tracer = pysimple.simple.Tracer()
        
        # Initialize field
        pysimple.simple_main.init_field(tracer, 'wout.nc', 5, 5, 5, 1)
        
        # Evaluate field at coordinates
        npoints = len(coords)
        results = {
            'modB': np.zeros(npoints),
            'sqrtg': np.zeros(npoints),
            'dBds': np.zeros(npoints),
            'dBdtheta': np.zeros(npoints),
            'dBdzeta': np.zeros(npoints),
            'Bs_cov': np.zeros(npoints),
            'Btheta_cov': np.zeros(npoints),
            'Bzeta_cov': np.zeros(npoints),
        }
        
        # Temporary arrays for field evaluation
        bder = np.empty(3)
        hcovar = np.empty(3)
        hctrvr = np.empty(3)
        hcurl = np.empty(3)
        
        for i, (s, theta, zeta) in enumerate(coords):
            x = np.array([s, theta, zeta])
            
            try:
                # Evaluate field using magfie routine
                if field_type == 2:  # Boozer
                    modB, sqrtg = pysimple.magfie_sub.magfie_boozer(x, bder, hcovar, hctrvr, hcurl)
                elif field_type == 5:  # BOOZXFORM
                    modB, sqrtg = pysimple.magfie_sub.magfie_boozxform(x, bder, hcovar, hctrvr, hcurl)
                else:
                    raise ValueError(f"Unsupported field type: {field_type}")
                
                results['modB'][i] = modB
                results['sqrtg'][i] = sqrtg
                results['dBds'][i] = bder[0] * modB
                results['dBdtheta'][i] = bder[1] * modB
                results['dBdzeta'][i] = bder[2] * modB
                results['Bs_cov'][i] = hcovar[0] * modB
                results['Btheta_cov'][i] = hcovar[1] * modB
                results['Bzeta_cov'][i] = hcovar[2] * modB
                
            except Exception as e:
                print(f"Error evaluating field at {x}: {e}")
                # Fill with NaN on error
                for key in results:
                    results[key][i] = np.nan
        
        return results
        
    finally:
        # Clean up temporary file
        os.unlink(config_file)


def plot_field_comparison(coords, boozer_results, boozxform_results):
    """Create comparison plots of field evaluation results"""
    
    s_vals = [c[0] for c in coords]
    theta_vals = [c[1] for c in coords]
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle('SIMPLE Boozer vs BOOZXFORM Field Evaluation Comparison', fontsize=14)
    
    # Plot 1: |B| comparison
    ax = axes[0, 0]
    ax.plot(theta_vals, boozer_results['modB'], 'b-', label='Internal Boozer', linewidth=2)
    ax.plot(theta_vals, boozxform_results['modB'], 'r--', label='BOOZXFORM', linewidth=2)
    ax.set_xlabel('θ_Boozer')
    ax.set_ylabel('|B|')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title('Magnetic Field Strength')
    
    # Plot 2: sqrt(g) comparison
    ax = axes[0, 1]
    ax.plot(theta_vals, boozer_results['sqrtg'], 'b-', label='Internal Boozer', linewidth=2)
    ax.plot(theta_vals, boozxform_results['sqrtg'], 'r--', label='BOOZXFORM', linewidth=2)
    ax.set_xlabel('θ_Boozer')
    ax.set_ylabel('√g')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title('Jacobian')
    
    # Plot 3: B derivatives
    ax = axes[0, 2]
    ax.plot(theta_vals, boozer_results['dBds'], 'b-', label='∂B/∂s (Boozer)', linewidth=2)
    ax.plot(theta_vals, boozxform_results['dBds'], 'r--', label='∂B/∂s (BOOZXFORM)', linewidth=2)
    ax.set_xlabel('θ_Boozer')
    ax.set_ylabel('∂B/∂s')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title('Field Derivative')
    
    # Plot 4: B covariant components
    ax = axes[1, 0]
    ax.plot(theta_vals, boozer_results['Btheta_cov'], 'b-', label='B_θ (Boozer)', linewidth=2)
    ax.plot(theta_vals, boozxform_results['Btheta_cov'], 'r--', label='B_θ (BOOZXFORM)', linewidth=2)
    ax.set_xlabel('θ_Boozer')
    ax.set_ylabel('B_θ')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title('Poloidal Field Component')
    
    # Plot 5: Relative differences
    ax = axes[1, 1]
    modB_diff = np.abs(boozer_results['modB'] - boozxform_results['modB']) / np.abs(boozer_results['modB'])
    sqrtg_diff = np.abs(boozer_results['sqrtg'] - boozxform_results['sqrtg']) / np.abs(boozer_results['sqrtg'])
    
    ax.semilogy(theta_vals, modB_diff, 'g-', label='|B| rel. diff', linewidth=2)
    ax.semilogy(theta_vals, sqrtg_diff, 'm-', label='√g rel. diff', linewidth=2)
    ax.set_xlabel('θ_Boozer')
    ax.set_ylabel('Relative Difference')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title('Relative Differences')
    
    # Plot 6: Summary statistics
    ax = axes[1, 2]
    ax.axis('off')
    
    # Calculate statistics
    modB_mean_diff = np.nanmean(modB_diff)
    modB_max_diff = np.nanmax(modB_diff)
    sqrtg_mean_diff = np.nanmean(sqrtg_diff)
    sqrtg_max_diff = np.nanmax(sqrtg_diff)
    
    valid_boozer = np.sum(~np.isnan(boozer_results['modB']))
    valid_boozxform = np.sum(~np.isnan(boozxform_results['modB']))
    
    stats_text = f"""Comparison Statistics:
    
Valid Evaluations:
  Internal Boozer: {valid_boozer}/{len(coords)}
  BOOZXFORM: {valid_boozxform}/{len(coords)}

|B| Differences:
  Mean: {modB_mean_diff:.2e}
  Max:  {modB_max_diff:.2e}

√g Differences:
  Mean: {sqrtg_mean_diff:.2e}
  Max:  {sqrtg_max_diff:.2e}

|B| Range:
  Boozer: [{np.nanmin(boozer_results['modB']):.1f}, {np.nanmax(boozer_results['modB']):.1f}]
  BOOZXFORM: [{np.nanmin(boozxform_results['modB']):.1f}, {np.nanmax(boozxform_results['modB']):.1f}]
"""
    
    ax.text(0.1, 0.9, stats_text, transform=ax.transAxes, fontsize=10, 
            verticalalignment='top', fontfamily='monospace')
    
    plt.tight_layout()
    return fig


def main():
    """Main comparison routine"""
    print("=== SIMPLE Field Evaluation Comparison ===")
    print("Comparing internal Boozer vs BOOZXFORM field evaluation")
    
    # Check files exist
    if not os.path.exists('wout.nc'):
        print("Error: wout.nc not found. Run setup_booz_xform_test.sh first.")
        return
    
    if not os.path.exists('boozmn_LandremanPaul2021_QA_lowres.nc'):
        print("Error: booz_xform file not found. Run setup_booz_xform_test.sh first.")
        return
    
    # Define test coordinates in Boozer space
    # Test at s=0.5, various poloidal angles, zeta=0
    s_test = 0.5
    theta_test = np.linspace(0, 2*np.pi, 20)
    zeta_test = 0.0
    
    coords = [(s_test, theta, zeta_test) for theta in theta_test]
    
    print(f"\nEvaluating field at {len(coords)} points:")
    print(f"  s = {s_test}")
    print(f"  θ ∈ [0, 2π] ({len(theta_test)} points)")
    print(f"  ζ = {zeta_test}")
    
    # Evaluate with internal Boozer
    print("\n1. Evaluating with internal Boozer conversion (isw_field_type=2)...")
    try:
        boozer_results = evaluate_field_at_coordinates(2, coords)
        print(f"   Internal Boozer: {np.sum(~np.isnan(boozer_results['modB']))}/{len(coords)} successful evaluations")
    except Exception as e:
        print(f"   Error with internal Boozer: {e}")
        boozer_results = None
    
    # Evaluate with BOOZXFORM
    print("\n2. Evaluating with BOOZXFORM file (isw_field_type=5)...")
    try:
        boozxform_results = evaluate_field_at_coordinates(5, coords, 'boozmn_LandremanPaul2021_QA_lowres.nc')
        print(f"   BOOZXFORM: {np.sum(~np.isnan(boozxform_results['modB']))}/{len(coords)} successful evaluations")
    except Exception as e:
        print(f"   Error with BOOZXFORM: {e}")
        boozxform_results = None
    
    # Create comparison plots
    if boozer_results is not None and boozxform_results is not None:
        print("\n3. Creating comparison plots...")
        fig = plot_field_comparison(coords, boozer_results, boozxform_results)
        
        # Save figure
        output_file = 'field_evaluation_comparison.png'
        fig.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"   Plots saved to: {output_file}")
        
        # Print summary
        print("\n=== Comparison Summary ===")
        if np.any(~np.isnan(boozer_results['modB'])) and np.any(~np.isnan(boozxform_results['modB'])):
            valid_mask = ~np.isnan(boozer_results['modB']) & ~np.isnan(boozxform_results['modB'])
            if np.any(valid_mask):
                modB_rel_diff = np.abs(boozer_results['modB'][valid_mask] - boozxform_results['modB'][valid_mask]) / np.abs(boozer_results['modB'][valid_mask])
                print(f"Field strength relative difference: mean={np.mean(modB_rel_diff):.2e}, max={np.max(modB_rel_diff):.2e}")
            else:
                print("No valid comparisons available")
        else:
            print("Insufficient data for meaningful comparison")
    else:
        print("\n3. Cannot create comparison plots - field evaluation failed")
        if boozer_results is None:
            print("   Internal Boozer evaluation failed")
        if boozxform_results is None:
            print("   BOOZXFORM evaluation failed")


if __name__ == '__main__':
    main()