#!/usr/bin/env python3
"""
Plot small-s behavior of VMEC and GVEC fields
"""

import numpy as np
import matplotlib.pyplot as plt

import os

# Load the data from build directory
data_dir = os.path.join(os.path.dirname(__file__), '../../build/test/tests')
vmec_data = np.loadtxt(os.path.join(data_dir, 'Bmod_vmec_radial.dat'), skiprows=1)
gvec_data = np.loadtxt(os.path.join(data_dir, 'Bmod_gvec_radial.dat'), skiprows=1)

# Extract columns
s_vmec = vmec_data[:, 0]
B_vmec = vmec_data[:, 1]
s_gvec = gvec_data[:, 0]
B_gvec = gvec_data[:, 1]

# Focus on small s region
small_s_mask_vmec = s_vmec < 0.2
small_s_mask_gvec = s_gvec < 0.2

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle('Small-s Behavior (s < 0.2)', fontsize=14)

# Linear plot
ax1.plot(s_vmec[small_s_mask_vmec], B_vmec[small_s_mask_vmec], 'b-', linewidth=2, label='VMEC')
ax1.plot(s_gvec[small_s_mask_gvec], B_gvec[small_s_mask_gvec], 'r--', linewidth=2, label='GVEC')
ax1.set_xlabel('Normalized flux s')
ax1.set_ylabel('Magnetic field magnitude |B| [T]')
ax1.grid(True, alpha=0.3)
ax1.legend()

# Log-log plot
mask_vmec = (s_vmec > 1e-3) & (s_vmec < 0.2) & (B_vmec > 0)
mask_gvec = (s_gvec > 1e-3) & (s_gvec < 0.2) & (B_gvec > 0)

ax2.loglog(s_vmec[mask_vmec], B_vmec[mask_vmec], 'b-', linewidth=2, label='VMEC')
ax2.loglog(s_gvec[mask_gvec], B_gvec[mask_gvec], 'r--', linewidth=2, label='GVEC')
ax2.set_xlabel('Normalized flux s')
ax2.set_ylabel('Magnetic field magnitude |B| [T]')
ax2.set_title('Small-s Behavior (log-log)')
ax2.grid(True, alpha=0.3, which='both')
ax2.legend()

plt.tight_layout()
output_dir = os.path.join(os.path.dirname(__file__), '../../build/test/tests')
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(output_dir, 'small_s_behavior.png')
plt.savefig(output_file, dpi=150, bbox_inches='tight')
print(f"Saved small-s behavior plot to {output_file}")

# Print axis values
print("\nAxis field values:")
idx_vmec_0 = np.argmin(np.abs(s_vmec))
idx_gvec_0 = np.argmin(np.abs(s_gvec))
print(f"VMEC at s={s_vmec[idx_vmec_0]:.4f}: B = {B_vmec[idx_vmec_0]:.1f} T")
print(f"GVEC at s={s_gvec[idx_gvec_0]:.4f}: B = {B_gvec[idx_gvec_0]:.1f} T")

# Check scaling behavior near axis
if len(s_vmec[mask_vmec]) > 10 and len(s_gvec[mask_gvec]) > 10:
    # Fit power law B ~ s^n for small s
    from scipy.optimize import curve_fit
    
    def power_law(s, A, n):
        return A * s**n
    
    # Fit VMEC
    s_fit = s_vmec[mask_vmec][:20]  # First 20 points
    B_fit = B_vmec[mask_vmec][:20]
    try:
        popt_vmec, _ = curve_fit(power_law, s_fit, B_fit)
        print(f"\nVMEC scaling near axis: B ~ {popt_vmec[0]:.1f} * s^{popt_vmec[1]:.3f}")
    except:
        print("\nCould not fit VMEC power law")
    
    # Fit GVEC
    s_fit = s_gvec[mask_gvec][:20]
    B_fit = B_gvec[mask_gvec][:20]
    try:
        popt_gvec, _ = curve_fit(power_law, s_fit, B_fit)
        print(f"GVEC scaling near axis: B ~ {popt_gvec[0]:.1f} * s^{popt_gvec[1]:.3f}")
    except:
        print("Could not fit GVEC power law")