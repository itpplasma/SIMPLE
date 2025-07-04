#!/usr/bin/env python3
"""
Create 1D radial comparison plots from the actual test data.
This script analyzes the real field comparison data.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def check_for_existing_data():
    """Check if test data files exist"""
    files_to_check = [
        'Bmod_vmec_2d.dat',
        'Bmod_gvec_2d.dat',
        'radial_comparison.dat'
    ]
    
    existing_files = []
    for file in files_to_check:
        if os.path.exists(file):
            existing_files.append(file)
            print(f"Found existing file: {file}")
    
    return existing_files

def extract_radial_profiles():
    """Extract radial profiles from 2D data"""
    print("Extracting radial profiles from 2D field data...")
    
    # Load VMEC data
    try:
        vmec_data = np.loadtxt('Bmod_vmec_2d.dat', comments='#')
        print(f"Loaded VMEC data: {vmec_data.shape}")
    except FileNotFoundError:
        print("VMEC data file not found")
        return None
    
    # Load GVEC data  
    try:
        gvec_data = np.loadtxt('Bmod_gvec_2d.dat', comments='#')
        print(f"Loaded GVEC data: {gvec_data.shape}")
    except FileNotFoundError:
        print("GVEC data file not found")
        return None
    
    # Extract columns: s, theta, R, Z, |B|
    s_vmec = vmec_data[:, 0]
    theta_vmec = vmec_data[:, 1]
    Bmod_vmec = vmec_data[:, 4]
    
    s_gvec = gvec_data[:, 0]
    theta_gvec = gvec_data[:, 1]
    Bmod_gvec = gvec_data[:, 4]
    
    # Extract radial profile at theta = 0
    # Find points with theta closest to 0
    theta_target = 0.0
    
    # For VMEC
    theta_diff_vmec = np.abs(theta_vmec - theta_target)
    vmec_mask = theta_diff_vmec < 0.1  # Within 0.1 radian of theta=0
    
    # For GVEC
    theta_diff_gvec = np.abs(theta_gvec - theta_target)
    gvec_mask = theta_diff_gvec < 0.1
    
    s_profile_vmec = s_vmec[vmec_mask]
    B_profile_vmec = Bmod_vmec[vmec_mask]
    
    s_profile_gvec = s_gvec[gvec_mask]
    B_profile_gvec = Bmod_gvec[gvec_mask]
    
    # Sort by s
    vmec_sort = np.argsort(s_profile_vmec)
    gvec_sort = np.argsort(s_profile_gvec)
    
    s_profile_vmec = s_profile_vmec[vmec_sort]
    B_profile_vmec = B_profile_vmec[vmec_sort]
    
    s_profile_gvec = s_profile_gvec[gvec_sort]
    B_profile_gvec = B_profile_gvec[gvec_sort]
    
    print(f"Extracted radial profiles:")
    print(f"  VMEC: {len(s_profile_vmec)} points")
    print(f"  GVEC: {len(s_profile_gvec)} points")
    
    return {
        's_vmec': s_profile_vmec,
        'B_vmec': B_profile_vmec,
        's_gvec': s_profile_gvec,
        'B_gvec': B_profile_gvec
    }

def interpolate_and_compare(profiles):
    """Interpolate profiles to common grid and compare"""
    print("Interpolating profiles to common grid...")
    
    # Create common s grid
    s_min = max(profiles['s_vmec'].min(), profiles['s_gvec'].min())
    s_max = min(profiles['s_vmec'].max(), profiles['s_gvec'].max())
    s_common = np.linspace(s_min, s_max, 100)
    
    # Interpolate both profiles
    B_vmec_interp = np.interp(s_common, profiles['s_vmec'], profiles['B_vmec'])
    B_gvec_interp = np.interp(s_common, profiles['s_gvec'], profiles['B_gvec'])
    
    # Calculate relative error
    rel_error = np.abs(B_gvec_interp - B_vmec_interp) / np.abs(B_vmec_interp)
    
    return {
        's_common': s_common,
        'B_vmec_interp': B_vmec_interp,
        'B_gvec_interp': B_gvec_interp,
        'rel_error': rel_error
    }

def create_1d_analysis_plots(profiles, comparison):
    """Create comprehensive 1D analysis plots"""
    print("Creating 1D analysis plots...")
    
    # Create figure with subplots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    
    # Plot 1: Raw profiles
    ax1.plot(profiles['s_vmec'], profiles['B_vmec'], 'bo-', label='VMEC', markersize=3)
    ax1.plot(profiles['s_gvec'], profiles['B_gvec'], 'ro-', label='GVEC', markersize=3)
    ax1.set_xlabel('Normalized flux s')
    ax1.set_ylabel('Magnetic field magnitude |B| [T]')
    ax1.set_title('Raw Field Profiles at θ ≈ 0')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Interpolated comparison
    ax2.plot(comparison['s_common'], comparison['B_vmec_interp'], 'b-', label='VMEC', linewidth=2)
    ax2.plot(comparison['s_common'], comparison['B_gvec_interp'], 'r--', label='GVEC', linewidth=2)
    ax2.set_xlabel('Normalized flux s')
    ax2.set_ylabel('Magnetic field magnitude |B| [T]')
    ax2.set_title('Interpolated Field Comparison')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Log-log comparison
    ax3.loglog(comparison['s_common'], comparison['B_vmec_interp'], 'b-', label='VMEC', linewidth=2)
    ax3.loglog(comparison['s_common'], comparison['B_gvec_interp'], 'r--', label='GVEC', linewidth=2)
    ax3.set_xlabel('Normalized flux s')
    ax3.set_ylabel('Magnetic field magnitude |B| [T]')
    ax3.set_title('Log-Log Field Comparison')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Relative error
    ax4.semilogy(comparison['s_common'], comparison['rel_error'], 'g-', linewidth=2)
    ax4.set_xlabel('Normalized flux s')
    ax4.set_ylabel('Relative error |B_gvec - B_vmec| / |B_vmec|')
    ax4.set_title('Relative Error vs Flux')
    ax4.grid(True, alpha=0.3)
    ax4.axhline(y=0.01, color='r', linestyle='--', alpha=0.7, label='1% error')
    ax4.axhline(y=0.1, color='orange', linestyle='--', alpha=0.7, label='10% error')
    ax4.legend()
    
    plt.tight_layout()
    plt.savefig('actual_1d_field_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create detailed small-s analysis
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Focus on small s
    small_s_mask = comparison['s_common'] < 0.2
    s_small = comparison['s_common'][small_s_mask]
    B_vmec_small = comparison['B_vmec_interp'][small_s_mask]
    B_gvec_small = comparison['B_gvec_interp'][small_s_mask]
    
    # Plot 1: Small-s linear
    ax1.plot(s_small, B_vmec_small, 'b-', label='VMEC', linewidth=2)
    ax1.plot(s_small, B_gvec_small, 'r--', label='GVEC', linewidth=2)
    ax1.set_xlabel('Normalized flux s')
    ax1.set_ylabel('Magnetic field magnitude |B| [T]')
    ax1.set_title('Small-s Behavior (s < 0.2)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Small-s log
    ax2.loglog(s_small, B_vmec_small, 'b-', label='VMEC', linewidth=2)
    ax2.loglog(s_small, B_gvec_small, 'r--', label='GVEC', linewidth=2)
    ax2.set_xlabel('Normalized flux s')
    ax2.set_ylabel('Magnetic field magnitude |B| [T]')
    ax2.set_title('Small-s Behavior (log-log)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('actual_small_s_behavior.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("1D analysis plots created:")
    print("  - actual_1d_field_comparison.png")
    print("  - actual_small_s_behavior.png")
    
    return True

def analyze_field_behavior(comparison):
    """Analyze field behavior and print statistics"""
    print("\n" + "="*60)
    print("FIELD BEHAVIOR ANALYSIS")
    print("="*60)
    
    s_common = comparison['s_common']
    B_vmec = comparison['B_vmec_interp']
    B_gvec = comparison['B_gvec_interp']
    rel_error = comparison['rel_error']
    
    print(f"Flux range: s = {s_common.min():.4f} to {s_common.max():.4f}")
    print(f"Number of points: {len(s_common)}")
    print()
    
    print(f"Field magnitude range:")
    print(f"  VMEC: {B_vmec.min():.3f} to {B_vmec.max():.3f} T")
    print(f"  GVEC: {B_gvec.min():.3f} to {B_gvec.max():.3f} T")
    print()
    
    print(f"Relative error statistics:")
    print(f"  Mean: {rel_error.mean():.3f}")
    print(f"  Median: {np.median(rel_error):.3f}")
    print(f"  Max: {rel_error.max():.3f}")
    print(f"  Min: {rel_error.min():.3f}")
    print()
    
    # Small-s analysis
    small_s_mask = s_common < 0.1
    if np.any(small_s_mask):
        print(f"Small-s behavior (s < 0.1):")
        print(f"  VMEC field range: {B_vmec[small_s_mask].min():.3f} to {B_vmec[small_s_mask].max():.3f} T")
        print(f"  GVEC field range: {B_gvec[small_s_mask].min():.3f} to {B_gvec[small_s_mask].max():.3f} T")
        print(f"  Average field ratio: {(B_gvec[small_s_mask] / B_vmec[small_s_mask]).mean():.3f}")
        print(f"  Max field ratio: {(B_gvec[small_s_mask] / B_vmec[small_s_mask]).max():.3f}")
    print()
    
    # Very small-s analysis
    very_small_s_mask = s_common < 0.01
    if np.any(very_small_s_mask):
        print(f"Very small-s behavior (s < 0.01):")
        print(f"  VMEC field: {B_vmec[very_small_s_mask].mean():.3f} T (average)")
        print(f"  GVEC field: {B_gvec[very_small_s_mask].mean():.3f} T (average)")
        print(f"  Ratio: {(B_gvec[very_small_s_mask] / B_vmec[very_small_s_mask]).mean():.3f}")
    print()
    
    # Identify worst regions
    worst_error_idx = np.argmax(rel_error)
    print(f"Worst error region:")
    print(f"  s = {s_common[worst_error_idx]:.4f}")
    print(f"  VMEC field: {B_vmec[worst_error_idx]:.3f} T")
    print(f"  GVEC field: {B_gvec[worst_error_idx]:.3f} T")
    print(f"  Relative error: {rel_error[worst_error_idx]:.3f}")

def main():
    """Main function"""
    print("1D RADIAL FIELD COMPARISON ANALYSIS")
    print("="*50)
    
    # Check for existing data
    existing_files = check_for_existing_data()
    
    if not existing_files:
        print("No test data files found. Please run the field comparison test first.")
        print("Expected files: Bmod_vmec_2d.dat, Bmod_gvec_2d.dat")
        return 1
    
    # Extract radial profiles
    profiles = extract_radial_profiles()
    if profiles is None:
        return 1
    
    # Interpolate and compare
    comparison = interpolate_and_compare(profiles)
    
    # Create analysis plots
    create_1d_analysis_plots(profiles, comparison)
    
    # Analyze behavior
    analyze_field_behavior(comparison)
    
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print("Successfully created 1D radial comparison plots from actual test data.")
    print("The plots show the detailed radial behavior of both field implementations.")
    print("Key findings should be visible in the small-s behavior plots.")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())