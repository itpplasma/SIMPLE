#!/usr/bin/env python3
"""
Compute the exact scaling factor between VMEC and GVEC fields
"""
import numpy as np

# Load the 2D data
vmec_data = np.loadtxt('build/test/tests/Bmod_vmec_2d.dat', comments='#')
gvec_data = np.loadtxt('build/test/tests/Bmod_gvec_2d.dat', comments='#')

# Extract field magnitudes
B_vmec = vmec_data[:, 4]
B_gvec = gvec_data[:, 4]

# Remove very small values to avoid division issues
mask = (B_vmec > 1e-10) & (B_gvec > 1e-10)
B_vmec_valid = B_vmec[mask]
B_gvec_valid = B_gvec[mask]

# Compute ratios
ratios = B_vmec_valid / B_gvec_valid

# Statistics
mean_ratio = np.mean(ratios)
median_ratio = np.median(ratios)
std_ratio = np.std(ratios)

print("VMEC/GVEC Field Scaling Analysis")
print("="*50)
print(f"Mean ratio: {mean_ratio:.6f}")
print(f"Median ratio: {median_ratio:.6f}")
print(f"Std deviation: {std_ratio:.6f}")
print(f"Min ratio: {np.min(ratios):.6f}")
print(f"Max ratio: {np.max(ratios):.6f}")
print()
print(f"Recommended scaling factor: {median_ratio:.1f}")
print()
print("To make GVEC match VMEC, multiply GVEC by this factor.")