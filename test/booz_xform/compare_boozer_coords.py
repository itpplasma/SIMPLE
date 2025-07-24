#!/usr/bin/env python3
"""
Compare Boozer coordinates from SIMPLE's internal conversion 
vs booz_xform external file
"""
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import sys
import os

# Add parent directory to path for pysimple
# Handle both source tree and build tree execution
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


def read_booz_xform_file(filename):
    """Read BOOZXFORM NetCDF file and return key data"""
    with nc.Dataset(filename, 'r') as f:
        # This is the older BOOZXFORM format from Stellopt
        data = {
            'ns': int(f.variables['ns_b'][:]),
            'nfp': int(f.variables['nfp_b'][:]),
            'mboz': int(f.variables['mboz_b'][:]),
            'nboz': int(f.variables['nboz_b'][:]),
            'mnboz': int(f.variables['mnboz_b'][:]),
            'iota': f.variables['iota_b'][:],
            'buco': f.variables['buco_b'][:],  # Boozer I (poloidal covariant)
            'bvco': f.variables['bvco_b'][:],  # Boozer G (toroidal covariant)
            'beta': f.variables['beta_b'][:],
            'phip': f.variables['phip_b'][:],
            'chi': f.variables['chi_b'][:],
            'pres': f.variables['pres_b'][:],
            'ixm': f.variables['ixm_b'][:].astype(int),
            'ixn': f.variables['ixn_b'][:].astype(int),
            'rmnc': f.variables['rmnc_b'][:],
            'zmns': f.variables['zmns_b'][:],
            'pmns': f.variables['pmns_b'][:],  # p = theta_VMEC - theta_Boozer
            'gmn': f.variables['gmn_b'][:],    # g = zeta_VMEC - zeta_Boozer
            'bmnc': f.variables['bmnc_b'][:],   # |B| Fourier coefficients
            'jlist': f.variables['jlist'][:].astype(int),
        }
        
        # Create s array (normalized toroidal flux)
        data['s'] = np.linspace(0, 1, data['ns'])
        
        # Check for stellarator symmetry
        try:
            data['lasym'] = bool(f.variables['lasym__logical__'][:])
        except:
            data['lasym'] = False
            
        print(f"Loaded BOOZXFORM file: {filename}")
        print(f"  ns = {data['ns']}, nfp = {data['nfp']}")
        print(f"  mboz = {data['mboz']}, nboz = {data['nboz']}")
        print(f"  mnboz = {data['mnboz']} Fourier modes")
        print(f"  stellarator symmetric: {not data['lasym']}")
        
    return data


def setup_simple_boozer(vmec_file):
    """Initialize SIMPLE with internal Boozer conversion"""
    tracer = pysimple.simple.Tracer()
    
    # Read config
    pysimple.params.read_config("simple.in")
    
    # Set to use internal Boozer (isw_field_type = 2)
    pysimple.velo_mod.isw_field_type = 2
    
    # Initialize field
    pysimple.simple_main.init_field(tracer, vmec_file, 5, 5, 5, 1)
    
    return tracer


def evaluate_simple_boozer(s_vals, theta_vals, zeta_vals):
    """Evaluate SIMPLE's internal Boozer coordinates at given points"""
    npoints = len(s_vals)
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
    
    for i in range(npoints):
        x = np.array([s_vals[i], theta_vals[i], zeta_vals[i]])
        
        # Evaluate field using SIMPLE's Boozer routine
        modB, sqrtg = pysimple.magfie_sub.magfie_boozer(x, bder, hcovar, hctrvr, hcurl)
        
        results['modB'][i] = modB
        results['sqrtg'][i] = sqrtg
        results['dBds'][i] = bder[0] * modB
        results['dBdtheta'][i] = bder[1] * modB
        results['dBdzeta'][i] = bder[2] * modB
        results['Bs_cov'][i] = hcovar[0] * modB
        results['Btheta_cov'][i] = hcovar[1] * modB
        results['Bzeta_cov'][i] = hcovar[2] * modB
        
    return results


def compare_coordinate_transform():
    """Compare VMEC->Boozer transformation from SIMPLE vs booz_xform"""
    
    # Test points in VMEC coordinates
    s_test = np.array([0.25, 0.5, 0.75])
    theta_vmec = np.linspace(0, 2*np.pi, 20)
    zeta_vmec = np.array([0.0])  # Single toroidal angle
    
    # Arrays to store results
    theta_boozer_simple = np.zeros((len(s_test), len(theta_vmec)))
    zeta_boozer_simple = np.zeros((len(s_test), len(theta_vmec)))
    
    # Check if boozer conversion is available
    try:
        # Try to access boozer conversion
        if hasattr(pysimple, 'boozer_sub'):
            boozer_sub = pysimple.boozer_sub
        elif hasattr(pysimple, 'Boozer_Coordinates_Mod'):
            # Try through Boozer_Coordinates_Mod
            boozer_sub = pysimple.Boozer_Coordinates_Mod
        elif hasattr(pysimple, 'boozer_coordinates_mod'):
            # Try lowercase version
            boozer_sub = pysimple.boozer_coordinates_mod
        else:
            print("Warning: Boozer conversion not directly accessible")
            print("Available modules:", [m for m in dir(pysimple) if not m.startswith('_')])
            # For now, just use the VMEC angles as a placeholder
            for i in range(len(s_test)):
                theta_boozer_simple[i, :] = theta_vmec
                zeta_boozer_simple[i, :] = zeta_vmec[0]
            return s_test, theta_vmec, theta_boozer_simple, zeta_boozer_simple
            
        # Get Boozer angles from SIMPLE's internal conversion
        for i, s in enumerate(s_test):
            for j, th_v in enumerate(theta_vmec):
                # Convert VMEC -> Boozer using SIMPLE
                if hasattr(boozer_sub, 'vmec_to_boozer'):
                    th_b, z_b = boozer_sub.vmec_to_boozer(s, th_v, zeta_vmec[0])
                else:
                    # Fallback - no conversion available
                    th_b, z_b = th_v, zeta_vmec[0]
                theta_boozer_simple[i, j] = th_b
                zeta_boozer_simple[i, j] = z_b
    except Exception as e:
        print(f"Error in coordinate conversion: {e}")
        # Use VMEC angles as fallback
        for i in range(len(s_test)):
            theta_boozer_simple[i, :] = theta_vmec
            zeta_boozer_simple[i, :] = zeta_vmec[0]
    
    return s_test, theta_vmec, theta_boozer_simple, zeta_boozer_simple


def plot_comparison(booz_data, s_test, theta_vmec, theta_boozer_simple):
    """Create comparison plots"""
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('SIMPLE Internal Boozer vs BOOZXFORM Comparison', fontsize=14)
    
    # Plot 1: Iota profile
    ax = axes[0, 0]
    ax.plot(booz_data['s'], booz_data['iota'], 'b-', label='BOOZXFORM', linewidth=2)
    # TODO: Add SIMPLE's iota when accessible
    ax.set_xlabel('s')
    ax.set_ylabel('ι (rotational transform)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title('Rotational Transform Profile')
    
    # Plot 2: Boozer coordinate transformation
    ax = axes[0, 1]
    colors = ['red', 'green', 'blue']
    for i, (s, color) in enumerate(zip(s_test, colors)):
        p_simple = theta_vmec - theta_boozer_simple[i, :]
        ax.plot(theta_vmec, p_simple, '-', color=color,
                label=f's={s:.2f} (SIMPLE)', linewidth=2)
    
    # Add BOOZXFORM p for comparison at s=0.5
    # We need to reconstruct p from Fourier coefficients
    s_idx = np.argmin(np.abs(booz_data['s'] - 0.5))
    if s_idx > 0 and s_idx < len(booz_data['jlist']):
        j_idx = booz_data['jlist'][s_idx-1]  # jlist is 1-indexed
        theta_test = theta_vmec
        p_booz = np.zeros_like(theta_test)
        # Sum Fourier components
        for k in range(booz_data['mnboz']):
            m = booz_data['ixm'][k]
            n = booz_data['ixn'][k]
            if j_idx >= 0 and j_idx < booz_data['pmns'].shape[0]:
                p_booz += booz_data['pmns'][j_idx, k] * np.sin(m * theta_test)
        ax.plot(theta_vmec, p_booz, 'k--', label='s=0.5 (BOOZXFORM)', linewidth=2)
    
    ax.set_xlabel('θ_VMEC')
    ax.set_ylabel('p = θ_VMEC - θ_Boozer')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title('Poloidal Angle Transform (ζ = 0)')
    
    # Plot 3: Fourier spectrum of |B|
    ax = axes[1, 0]
    # Show first few Fourier modes
    modes_to_show = 20
    mode_labels = [f'({booz_data["ixm"][k]},{booz_data["ixn"][k]})' 
                   for k in range(min(modes_to_show, booz_data['mnboz']))]
    s_idx = len(booz_data['jlist']) // 2
    if s_idx < booz_data['bmnc'].shape[0]:
        mode_amplitudes = np.abs(booz_data['bmnc'][s_idx, :modes_to_show])
        ax.bar(range(len(mode_amplitudes)), mode_amplitudes)
        ax.set_xlabel('Mode index')
        ax.set_ylabel('|B| Fourier amplitude')
        ax.set_title(f'|B| Fourier Spectrum at s≈{booz_data["s"][s_idx+1]:.2f}')
        ax.set_yscale('log')
    
    # Plot 4: G and I profiles
    ax = axes[1, 1]
    ax.plot(booz_data['s'], booz_data['bvco'], 'b-', label='G (toroidal)', linewidth=2)
    ax.plot(booz_data['s'], booz_data['buco'], 'r-', label='I (poloidal)', linewidth=2)
    ax.set_xlabel('s')
    ax.set_ylabel('Boozer covariant components')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title('Boozer I and G')
    
    plt.tight_layout()
    return fig


def main():
    """Main comparison routine"""
    print("=== Boozer Coordinate Comparison ===")
    
    # Check files exist
    if not os.path.exists('wout.nc'):
        print("Error: wout.nc not found. Run setup_booz_xform_test.sh first.")
        return
    
    if not os.path.exists('boozmn_LandremanPaul2021_QA_lowres.nc'):
        print("Error: booz_xform file not found. Run setup_booz_xform_test.sh first.")
        return
    
    # Load booz_xform data
    print("\n1. Loading booz_xform data...")
    booz_data = read_booz_xform_file('boozmn_LandremanPaul2021_QA_lowres.nc')
    
    # Initialize SIMPLE
    print("\n2. Initializing SIMPLE with internal Boozer conversion...")
    tracer = setup_simple_boozer('wout.nc')
    
    # Compare coordinate transformations
    print("\n3. Computing coordinate transformations...")
    s_test, theta_vmec, theta_boozer_simple, zeta_boozer_simple = compare_coordinate_transform()
    
    # Evaluate field quantities
    print("\n4. Evaluating field quantities...")
    # Test at mid-radius, various poloidal angles
    s_eval = 0.5 * np.ones(20)
    theta_eval = np.linspace(0, 2*np.pi, 20)
    zeta_eval = np.zeros(20)
    
    simple_results = evaluate_simple_boozer(s_eval, theta_eval, zeta_eval)
    
    # Create comparison plots
    print("\n5. Creating comparison plots...")
    fig = plot_comparison(booz_data, s_test, theta_vmec, theta_boozer_simple)
    
    # Save figure
    output_file = 'boozer_comparison.png'
    fig.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\nPlots saved to: {output_file}")
    
    # Print some numerical comparisons
    print("\n=== Numerical Comparison at s=0.5 ===")
    print(f"SIMPLE |B| range: [{simple_results['modB'].min():.4f}, {simple_results['modB'].max():.4f}]")
    print(f"SIMPLE sqrt(g): {simple_results['sqrtg'][0]:.6e}")
    
    # plt.show()  # Don't show in test mode


if __name__ == '__main__':
    main()