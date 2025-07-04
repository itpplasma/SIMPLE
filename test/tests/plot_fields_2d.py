#!/usr/bin/env python3
"""
Plot 2D magnetic field magnitude comparison between VMEC and GVEC
Reads data exported by export_field_2d.f90 and creates contour plots
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
import os
import sys

def read_grid_info(filename='grid_info.dat'):
    """Read grid information from Fortran output"""
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
            
        # Parse grid dimensions
        ns, nt = map(int, lines[0].split())
        
        # Parse grid ranges
        s_min, s_max, theta_min, theta_max = map(float, lines[1].split())
        
        # Parse fixed phi
        phi_fixed = float(lines[2].split()[0])
        
        return ns, nt, s_min, s_max, theta_min, theta_max, phi_fixed
        
    except FileNotFoundError:
        print(f"Error: {filename} not found. Run export_field_2d.x first.")
        sys.exit(1)

def read_field_data(filename):
    """Read 2D field data from Fortran output"""
    try:
        data = np.loadtxt(filename, comments='#')
        return data
    except FileNotFoundError:
        print(f"Error: {filename} not found. Run export_field_2d.x first.")
        sys.exit(1)

def create_2d_grid(data, ns, nt):
    """Reshape 1D data array to 2D grid"""
    # Remove empty lines that separate blocks in the output
    # Each block has nt points, followed by empty line
    valid_data = data[~np.isnan(data).any(axis=1)]
    
    # Reshape to 2D grid
    s_vals = valid_data[:, 0].reshape(ns, nt)
    theta_vals = valid_data[:, 1].reshape(ns, nt)
    R_vals = valid_data[:, 2].reshape(ns, nt)
    Z_vals = valid_data[:, 3].reshape(ns, nt)
    B_vals = valid_data[:, 4].reshape(ns, nt)
    
    return s_vals, theta_vals, R_vals, Z_vals, B_vals

def plot_field_comparison():
    """Create side-by-side contour plots of VMEC and GVEC fields"""
    
    print("Reading grid information...")
    ns, nt, s_min, s_max, theta_min, theta_max, phi_fixed = read_grid_info()
    
    print(f"Grid: {ns} x {nt} points")
    print(f"s range: {s_min:.3f} to {s_max:.3f}")
    print(f"theta range: {theta_min:.3f} to {theta_max:.3f}")
    print(f"phi fixed at: {phi_fixed:.3f}")
    
    print("\nReading VMEC field data...")
    vmec_data = read_field_data('Bmod_vmec_2d.dat')
    
    print("Reading GVEC field data...")
    gvec_data = read_field_data('Bmod_gvec_2d.dat')
    
    # Reshape data to 2D grids
    print("Reshaping data to 2D grids...")
    s_vmec, theta_vmec, R_vmec, Z_vmec, B_vmec = create_2d_grid(vmec_data, ns, nt)
    s_gvec, theta_gvec, R_gvec, Z_gvec, B_gvec = create_2d_grid(gvec_data, ns, nt)
    
    # Print field statistics
    print(f"\nField statistics:")
    print(f"VMEC |B| range: {B_vmec.min():.2e} to {B_vmec.max():.2e}")
    print(f"GVEC |B| range: {B_gvec.min():.2e} to {B_gvec.max():.2e}")
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle(f'Magnetic Field Magnitude Comparison (φ = {phi_fixed/np.pi:.2f}π)', fontsize=16)
    
    # Determine common color scale for better comparison
    vmin = min(B_vmec.min(), B_gvec.min())
    vmax = max(B_vmec.max(), B_gvec.max())
    
    # Use logarithmic scale if the range is large
    if vmax / vmin > 100:
        norm = colors.LogNorm(vmin=vmin, vmax=vmax)
        print("Using logarithmic color scale")
    else:
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
        print("Using linear color scale")
    
    # Plot 1: VMEC in (s, theta) coordinates
    ax1 = axes[0, 0]
    im1 = ax1.contourf(s_vmec, theta_vmec, B_vmec, levels=50, norm=norm, cmap='plasma')
    ax1.set_xlabel('s (normalized flux)')
    ax1.set_ylabel('θ (rad)')
    ax1.set_title('VMEC |B| in (s, θ)')
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: GVEC in (s, theta) coordinates
    ax2 = axes[0, 1]
    im2 = ax2.contourf(s_gvec, theta_gvec, B_gvec, levels=50, norm=norm, cmap='plasma')
    ax2.set_xlabel('s (normalized flux)')
    ax2.set_ylabel('θ (rad)')
    ax2.set_title('GVEC |B| in (s, θ)')
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: VMEC in (R, Z) coordinates (approximate)
    ax3 = axes[1, 0]
    im3 = ax3.contourf(R_vmec, Z_vmec, B_vmec, levels=50, norm=norm, cmap='plasma')
    ax3.set_xlabel('R (m)')
    ax3.set_ylabel('Z (m)')
    ax3.set_title('VMEC |B| in (R, Z) [approx]')
    ax3.set_aspect('equal')
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: GVEC in (R, Z) coordinates (approximate)
    ax4 = axes[1, 1]
    im4 = ax4.contourf(R_gvec, Z_gvec, B_gvec, levels=50, norm=norm, cmap='plasma')
    ax4.set_xlabel('R (m)')
    ax4.set_ylabel('Z (m)')
    ax4.set_title('GVEC |B| in (R, Z) [approx]')
    ax4.set_aspect('equal')
    ax4.grid(True, alpha=0.3)
    
    # Add colorbar
    cbar = fig.colorbar(im1, ax=axes, orientation='horizontal', 
                       fraction=0.05, pad=0.1, aspect=40)
    cbar.set_label('|B| (field units)', fontsize=12)
    
    plt.tight_layout()
    plt.savefig('field_comparison_2d.png', dpi=300, bbox_inches='tight')
    print(f"\nSaved comparison plot to: field_comparison_2d.png")
    
    # Create difference plot
    fig2, ax = plt.subplots(1, 1, figsize=(10, 8))
    
    # Calculate relative difference
    rel_diff = np.abs(B_gvec - B_vmec) / B_vmec * 100  # Percentage
    
    im = ax.contourf(s_vmec, theta_vmec, rel_diff, levels=50, cmap='RdYlBu_r')
    ax.set_xlabel('s (normalized flux)')
    ax.set_ylabel('θ (rad)')
    ax.set_title('Relative Difference: |B_GVEC - B_VMEC| / B_VMEC (%)')
    ax.grid(True, alpha=0.3)
    
    cbar2 = plt.colorbar(im, ax=ax)
    cbar2.set_label('Relative difference (%)', fontsize=12)
    
    plt.tight_layout()
    plt.savefig('field_difference_2d.png', dpi=300, bbox_inches='tight')
    print(f"Saved difference plot to: field_difference_2d.png")
    
    # Print summary statistics
    print(f"\nSummary statistics:")
    print(f"Maximum relative difference: {rel_diff.max():.2f}%")
    print(f"Mean relative difference: {rel_diff.mean():.2f}%")
    print(f"RMS relative difference: {np.sqrt(np.mean(rel_diff**2)):.2f}%")
    
    return B_vmec, B_gvec, rel_diff

def plot_individual_fields():
    """Create individual high-quality plots for each field"""
    
    print("\nCreating individual field plots...")
    
    # Read data
    ns, nt, s_min, s_max, theta_min, theta_max, phi_fixed = read_grid_info()
    vmec_data = read_field_data('Bmod_vmec_2d.dat')
    gvec_data = read_field_data('Bmod_gvec_2d.dat')
    
    # Reshape data
    s_vmec, theta_vmec, R_vmec, Z_vmec, B_vmec = create_2d_grid(vmec_data, ns, nt)
    s_gvec, theta_gvec, R_gvec, Z_gvec, B_gvec = create_2d_grid(gvec_data, ns, nt)
    
    # VMEC plot
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    im = ax.contourf(s_vmec, theta_vmec, B_vmec, levels=50, cmap='plasma')
    ax.set_xlabel('s (normalized flux)', fontsize=14)
    ax.set_ylabel('θ (rad)', fontsize=14)
    ax.set_title(f'VMEC Magnetic Field Magnitude (φ = {phi_fixed/np.pi:.2f}π)', fontsize=16)
    ax.grid(True, alpha=0.3)
    
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('|B| (VMEC internal units)', fontsize=12)
    
    plt.tight_layout()
    plt.savefig('vmec_field_2d.png', dpi=300, bbox_inches='tight')
    print("Saved VMEC field plot to: vmec_field_2d.png")
    plt.close()
    
    # GVEC plot
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    im = ax.contourf(s_gvec, theta_gvec, B_gvec, levels=50, cmap='plasma')
    ax.set_xlabel('s (normalized flux)', fontsize=14)
    ax.set_ylabel('θ (rad)', fontsize=14)
    ax.set_title(f'GVEC Magnetic Field Magnitude (φ = {phi_fixed/np.pi:.2f}π)', fontsize=16)
    ax.grid(True, alpha=0.3)
    
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('|B| (converted to VMEC units)', fontsize=12)
    
    plt.tight_layout()
    plt.savefig('gvec_field_2d.png', dpi=300, bbox_inches='tight')
    print("Saved GVEC field plot to: gvec_field_2d.png")
    plt.close()

if __name__ == "__main__":
    print("2D Magnetic Field Visualization")
    print("=" * 40)
    
    # Check if data files exist
    required_files = ['grid_info.dat', 'Bmod_vmec_2d.dat', 'Bmod_gvec_2d.dat']
    missing_files = [f for f in required_files if not os.path.exists(f)]
    
    if missing_files:
        print("Error: Missing data files:")
        for f in missing_files:
            print(f"  - {f}")
        print("\nPlease run export_field_2d.x first to generate the data files.")
        sys.exit(1)
    
    try:
        # Create comparison plots
        B_vmec, B_gvec, rel_diff = plot_field_comparison()
        
        # Create individual field plots
        plot_individual_fields()
        
        print("\n" + "=" * 40)
        print("Visualization completed successfully!")
        print("Generated files:")
        print("  - field_comparison_2d.png (side-by-side comparison)")
        print("  - field_difference_2d.png (relative difference)")
        print("  - vmec_field_2d.png (VMEC field only)")
        print("  - gvec_field_2d.png (GVEC field only)")
        
    except Exception as e:
        print(f"Error during visualization: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)