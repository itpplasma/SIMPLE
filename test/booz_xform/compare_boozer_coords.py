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
sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))

try:
    import pysimple
except ImportError:
    print("Error: pysimple not found. Please build SIMPLE first.")
    sys.exit(1)


def read_booz_xform_file(filename):
    """Read booz_xform NetCDF file and return key data"""
    with nc.Dataset(filename, 'r') as f:
        data = {
            'ns': f.dimensions['ns'].size,
            'nfp': int(f.variables['nfp'][:]),
            's': f.variables['s_b'][:],
            'iota': f.variables['iota'][:],
            'G': f.variables['G_b'][:],  # Boozer G (toroidal covariant component)
            'I': f.variables['I_b'][:],  # Boozer I (poloidal covariant component)
            'B': f.variables['B_b'][:],  # |B| on Boozer grid
            'dBds': f.variables['dBds_b'][:],
            'dBdtheta': f.variables['dBdtheta_b'][:],
            'dBdzeta': f.variables['dBdzeta_b'][:],
            'mpol': int(f.variables['mpol_b'][:]),
            'ntor': int(f.variables['ntor_b'][:]),
            'xm_b': f.variables['xm_b'][:],
            'xn_b': f.variables['xn_b'][:],
            'rmnc_b': f.variables['rmnc_b'][:],
            'zmns_b': f.variables['zmns_b'][:],
            'pmns_b': f.variables['pmns_b'][:],  # p = theta_VMEC - theta_Boozer
            'gmns_b': f.variables['gmns_b'][:],  # g = zeta_VMEC - zeta_Boozer
            'bmnc_b': f.variables['bmnc_b'][:],  # |B| Fourier coefficients
        }
        
        # Check for stellarator symmetry
        try:
            data['lasym'] = bool(f.variables['lasym__logical__'][:])
        except:
            data['lasym'] = False
            
        print(f"Loaded booz_xform file: {filename}")
        print(f"  ns = {data['ns']}, nfp = {data['nfp']}")
        print(f"  mpol = {data['mpol']}, ntor = {data['ntor']}")
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
    
    # Get Boozer angles from SIMPLE's internal conversion
    for i, s in enumerate(s_test):
        for j, th_v in enumerate(theta_vmec):
            # Convert VMEC -> Boozer using SIMPLE
            th_b, z_b = pysimple.boozer_sub.vmec_to_boozer(s, th_v, zeta_vmec[0])
            theta_boozer_simple[i, j] = th_b
            zeta_boozer_simple[i, j] = z_b
    
    # For comparison with booz_xform, we need to evaluate p and g
    # p = theta_VMEC - theta_Boozer, g = zeta_VMEC - zeta_Boozer
    
    return s_test, theta_vmec, theta_boozer_simple, zeta_boozer_simple


def plot_comparison(booz_data, s_test, theta_vmec, theta_boozer_simple):
    """Create comparison plots"""
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('SIMPLE Internal Boozer vs booz_xform Comparison', fontsize=14)
    
    # Plot 1: Iota profile
    ax = axes[0, 0]
    ax.plot(booz_data['s'], booz_data['iota'], 'b-', label='booz_xform', linewidth=2)
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
        ax.plot(theta_vmec, p_simple, f'{color}-', 
                label=f's={s:.2f} (SIMPLE)', linewidth=2)
    ax.set_xlabel('θ_VMEC')
    ax.set_ylabel('p = θ_VMEC - θ_Boozer')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title('Poloidal Angle Transform (ζ = 0)')
    
    # Plot 3: |B| on a surface
    ax = axes[1, 0]
    s_idx = len(booz_data['s']) // 2  # Middle surface
    if 'B' in booz_data and booz_data['B'].ndim > 1:
        # Plot |B| contours if available
        ntheta = booz_data['B'].shape[1] if booz_data['B'].ndim > 1 else 1
        nzeta = booz_data['B'].shape[2] if booz_data['B'].ndim > 2 else 1
        if ntheta > 1 and nzeta > 1:
            theta_grid = np.linspace(0, 2*np.pi, ntheta)
            zeta_grid = np.linspace(0, 2*np.pi/booz_data['nfp'], nzeta)
            THETA, ZETA = np.meshgrid(theta_grid, zeta_grid)
            cs = ax.contour(THETA, ZETA, booz_data['B'][s_idx, :, :].T)
            ax.clabel(cs, inline=True, fontsize=8)
            ax.set_title(f'|B| Contours at s={booz_data["s"][s_idx]:.2f}')
    else:
        # Just plot text if no 2D data
        ax.text(0.5, 0.5, f'|B| data shape: {booz_data["B"].shape}', 
                ha='center', va='center', transform=ax.transAxes)
        ax.set_title('|B| Data Structure')
    ax.set_xlabel('θ_Boozer')
    ax.set_ylabel('ζ_Boozer')
    
    # Plot 4: G and I profiles
    ax = axes[1, 1]
    ax.plot(booz_data['s'], booz_data['G'], 'b-', label='G (toroidal)', linewidth=2)
    ax.plot(booz_data['s'], booz_data['I'], 'r-', label='I (poloidal)', linewidth=2)
    ax.set_xlabel('s')
    ax.set_ylabel('Boozer covariant components')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title('Boozer G and I')
    
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
    
    plt.show()


if __name__ == '__main__':
    main()