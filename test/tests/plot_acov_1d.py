#!/usr/bin/env python3
"""
Plot 1D radial comparison of vector potential components between VMEC and GVEC
"""

import numpy as np
import matplotlib.pyplot as plt

# Read radial comparison data for |B|
bfield_data = np.loadtxt('radial_comparison.dat')
s = bfield_data[:, 0]
r = bfield_data[:, 1]
B_vmec = bfield_data[:, 2]
B_gvec = bfield_data[:, 3]

# Read vector potential data
try:
    acov_data = np.loadtxt('acov_comparison.dat')
    s_acov = acov_data[:, 0]
    r_acov = acov_data[:, 1]
    Acov1_vmec = acov_data[:, 2]
    Acov1_gvec = acov_data[:, 3]
    Acov2_vmec = acov_data[:, 4]
    Acov2_gvec = acov_data[:, 5]
    Acov3_vmec = acov_data[:, 6]
    Acov3_gvec = acov_data[:, 7]
    has_acov_data = True
except:
    print("Warning: acov_comparison.dat not found. Run test_vmec_gvec.x first.")
    has_acov_data = False

# Create figure with subplots
fig, axes = plt.subplots(4, 1, figsize=(10, 14))

# Plot 1: |B| comparison
ax1 = axes[0]
ax1.loglog(s, B_vmec, 'b-', label='VMEC', linewidth=2)
ax1.loglog(s, B_gvec, 'r--', label='GVEC', linewidth=2)
ax1.set_ylabel('|B| [Gauss]')
ax1.set_title('Magnetic Field Magnitude Comparison')
ax1.legend()
ax1.grid(True, alpha=0.3)

if has_acov_data:
    # Plot 2: Acov(1) - radial component
    ax2 = axes[1]
    mask_vmec = np.abs(Acov1_vmec) > 1e-10
    mask_gvec = np.abs(Acov1_gvec) > 1e-10
    
    if np.any(mask_vmec):
        ax2.semilogx(s_acov[mask_vmec], Acov1_vmec[mask_vmec]/1e8, 'b-', label='VMEC', linewidth=2)
    if np.any(mask_gvec):
        ax2.semilogx(s_acov[mask_gvec], Acov1_gvec[mask_gvec]/1e8, 'r--', label='GVEC', linewidth=2)
    
    ax2.set_ylabel('$A_s$ [$10^8$ Gauss·cm]')
    ax2.set_title('Vector Potential - Radial Component')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Acov(2) - poloidal component (main component)
    ax3 = axes[2]
    ax3.semilogx(s_acov, Acov2_vmec/1e8, 'b-', label='VMEC', linewidth=2)
    ax3.semilogx(s_acov, Acov2_gvec/1e8, 'r--', label='GVEC', linewidth=2)
    ax3.set_ylabel('$A_\\theta$ [$10^8$ Gauss·cm]')
    ax3.set_title('Vector Potential - Poloidal Component')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Acov(3) - toroidal component
    ax4 = axes[3]
    mask_vmec = np.abs(Acov3_vmec) > 1e-10
    mask_gvec = np.abs(Acov3_gvec) > 1e-10
    
    if np.any(mask_vmec):
        ax4.semilogx(s_acov[mask_vmec], Acov3_vmec[mask_vmec]/1e8, 'b-', label='VMEC', linewidth=2)
    if np.any(mask_gvec):
        ax4.semilogx(s_acov[mask_gvec], Acov3_gvec[mask_gvec]/1e8, 'r--', label='GVEC', linewidth=2)
    
    ax4.set_xlabel('s')
    ax4.set_ylabel('$A_\\phi$ [$10^8$ Gauss·cm]')
    ax4.set_title('Vector Potential - Toroidal Component')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
else:
    # Placeholder plots
    for i, (title, ylabel) in enumerate([
        ('Vector Potential - Radial Component', '$A_s$ [$10^8$ Gauss·cm]'),
        ('Vector Potential - Poloidal Component', '$A_\\theta$ [$10^8$ Gauss·cm]'),
        ('Vector Potential - Toroidal Component', '$A_\\phi$ [$10^8$ Gauss·cm]')
    ]):
        ax = axes[i+1]
        ax.text(0.5, 0.5, 'Run test_vmec_gvec.x to generate data', 
                transform=ax.transAxes, ha='center', va='center', fontsize=12)
        ax.set_xlabel('s')
        ax.set_ylabel(ylabel)
        ax.set_title(title)

plt.tight_layout()
plt.savefig('acov_comparison_1d.png', dpi=150, bbox_inches='tight')
print("Saved: acov_comparison_1d.png")

# Create a focused plot on Acov(2) comparison
if has_acov_data:
    fig2, ax = plt.subplots(figsize=(10, 6))
    
    # Linear plot to see the s-dependence better
    ax.plot(s_acov, Acov2_vmec/1e8, 'b-', label='VMEC', linewidth=2)
    ax.plot(s_acov, Acov2_gvec/1e8, 'r--', label='GVEC', linewidth=2)
    
    ax.set_xlabel('s')
    ax.set_ylabel('$A_\\theta$ [$10^8$ Gauss·cm]')
    ax.set_title('Vector Potential $A_\\theta$ - Linear Scale')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Add relative error as inset
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    axins = inset_axes(ax, width="40%", height="40%", loc='lower right')
    
    rel_error = np.abs(Acov2_gvec - Acov2_vmec) / np.abs(Acov2_vmec)
    axins.semilogy(s_acov, rel_error, 'g-', linewidth=2)
    axins.set_xlabel('s')
    axins.set_ylabel('Relative Error')
    axins.grid(True, alpha=0.3)
    
    plt.savefig('acov2_detailed_comparison.png', dpi=150, bbox_inches='tight')
    print("Saved: acov2_detailed_comparison.png")