#!/usr/bin/env python3
"""
Plot 1D radial comparison of VMEC and GVEC magnetic fields
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

# Create figure with 4 subplots
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle('VMEC vs GVEC Field Comparison', fontsize=14)

# 1. Raw field profiles
ax = axes[0, 0]
ax.plot(s_vmec, B_vmec, 'b-', linewidth=2, label='VMEC')
ax.plot(s_gvec, B_gvec, 'r--', linewidth=2, label='GVEC')
ax.set_xlabel('Normalized flux s')
ax.set_ylabel('Magnetic field magnitude |B| [T]')
ax.set_title('Raw Field Profiles at θ = 0')
ax.grid(True, alpha=0.3)
ax.legend()

# 2. Interpolated comparison
from scipy.interpolate import interp1d
f_gvec = interp1d(s_gvec, B_gvec, kind='cubic', bounds_error=False, fill_value='extrapolate')
B_gvec_interp = f_gvec(s_vmec)

ax = axes[0, 1]
ax.plot(s_vmec, B_vmec, 'b-', linewidth=2, label='VMEC')
ax.plot(s_vmec, B_gvec_interp, 'r--', linewidth=2, label='GVEC')
ax.set_xlabel('Normalized flux s')
ax.set_ylabel('Magnetic field magnitude |B| [T]')
ax.set_title('Interpolated Field Comparison')
ax.grid(True, alpha=0.3)
ax.legend()

# 3. Log-log plot
ax = axes[1, 0]
mask = (s_vmec > 1e-3) & (B_vmec > 0) & (B_gvec_interp > 0)
ax.loglog(s_vmec[mask], B_vmec[mask], 'b-', linewidth=2, label='VMEC')
ax.loglog(s_vmec[mask], B_gvec_interp[mask], 'r--', linewidth=2, label='GVEC')
ax.set_xlabel('Normalized flux s')
ax.set_ylabel('Magnetic field magnitude |B| [T]')
ax.set_title('Log-Log Field Comparison')
ax.grid(True, alpha=0.3, which='both')
ax.legend()

# 4. Relative error
ax = axes[1, 1]
rel_error = np.abs(B_gvec_interp - B_vmec) / B_vmec
ax.semilogy(s_vmec, rel_error, 'g-', linewidth=2, label='|B_gvec - B_vmec| / B_vmec')
ax.axhline(0.01, color='orange', linestyle='--', label='1% error')
ax.axhline(0.1, color='red', linestyle='--', label='10% error')
ax.set_xlabel('Normalized flux s')
ax.set_ylabel('Relative error')
ax.set_title('Relative Error vs Flux')
ax.grid(True, alpha=0.3)
ax.legend()
ax.set_ylim(1e-3, 1e0)

plt.tight_layout()
output_dir = os.path.join(os.path.dirname(__file__), '../../build/test/tests')
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(output_dir, 'field_comparison_1d.png')
plt.savefig(output_file, dpi=150, bbox_inches='tight')
print(f"Saved 1D field comparison plot to {output_file}")

# Print statistics
print("\nField comparison statistics:")
print(f"Mean relative error: {np.mean(rel_error):.3%}")
print(f"Max relative error: {np.max(rel_error):.3%}")
print(f"Error at s=0.5: {rel_error[np.argmin(np.abs(s_vmec - 0.5))]:.3%}")