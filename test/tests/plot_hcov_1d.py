#!/usr/bin/env python3
"""
Plot 1D radial profiles of hcov (normalized covariant B field) components 
comparing VMEC and GVEC implementations.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Read the radial comparison data
data = np.loadtxt('radial_comparison.dat')
s = data[:, 0]
r = data[:, 1]
Bmod_vmec = data[:, 2]
Bmod_gvec = data[:, 3]

# Read the hcov comparison data
hcov_data = np.loadtxt('hcov_comparison.dat')
s_h = hcov_data[:, 0]
r_h = hcov_data[:, 1]
hcov1_vmec = hcov_data[:, 2]
hcov1_gvec = hcov_data[:, 3]
hcov2_vmec = hcov_data[:, 4]
hcov2_gvec = hcov_data[:, 5]
hcov3_vmec = hcov_data[:, 6]
hcov3_gvec = hcov_data[:, 7]

# Create figure with subplots
fig = plt.figure(figsize=(12, 10))
gs = GridSpec(3, 2, figure=fig, hspace=0.3, wspace=0.3)

# Plot hcov(1) - radial component
ax1 = fig.add_subplot(gs[0, 0])
ax1.plot(s_h, hcov1_vmec, 'b-', label='VMEC', linewidth=2)
ax1.plot(s_h, hcov1_gvec, 'r--', label='GVEC', linewidth=2)
ax1.set_xlabel('s')
ax1.set_ylabel('hcov(1)')
ax1.set_title('Radial Component hcov(1)')
ax1.legend()
ax1.grid(True, alpha=0.3)
ax1.set_xlim([0, 1])

# Plot hcov(1) on log scale for small s
ax2 = fig.add_subplot(gs[0, 1])
mask = s_h > 1e-3
ax2.loglog(s_h[mask], np.abs(hcov1_vmec[mask]), 'b-', label='|hcov(1)| VMEC', linewidth=2)
ax2.loglog(s_h[mask], np.abs(hcov1_gvec[mask]), 'r--', label='|hcov(1)| GVEC', linewidth=2)
ax2.set_xlabel('s')
ax2.set_ylabel('|hcov(1)|')
ax2.set_title('Radial Component (log-log scale)')
ax2.legend()
ax2.grid(True, alpha=0.3, which='both')

# Plot hcov(2) - poloidal component
ax3 = fig.add_subplot(gs[1, 0])
ax3.plot(s_h, hcov2_vmec, 'b-', label='VMEC', linewidth=2)
ax3.plot(s_h, hcov2_gvec, 'r--', label='GVEC', linewidth=2)
ax3.set_xlabel('s')
ax3.set_ylabel('hcov(2)')
ax3.set_title('Poloidal Component hcov(2)')
ax3.legend()
ax3.grid(True, alpha=0.3)
ax3.set_xlim([0, 1])

# Plot hcov(2) relative error
ax4 = fig.add_subplot(gs[1, 1])
rel_error_2 = np.abs(hcov2_gvec - hcov2_vmec) / (np.abs(hcov2_vmec) + 1e-10)
ax4.semilogy(s_h, rel_error_2, 'g-', linewidth=2)
ax4.set_xlabel('s')
ax4.set_ylabel('Relative Error')
ax4.set_title('hcov(2) Relative Error')
ax4.grid(True, alpha=0.3)
ax4.set_xlim([0, 1])

# Plot hcov(3) - toroidal component
ax5 = fig.add_subplot(gs[2, 0])
ax5.plot(s_h, hcov3_vmec, 'b-', label='VMEC', linewidth=2)
ax5.plot(s_h, hcov3_gvec, 'r--', label='GVEC', linewidth=2)
ax5.set_xlabel('s')
ax5.set_ylabel('hcov(3)')
ax5.set_title('Toroidal Component hcov(3)')
ax5.legend()
ax5.grid(True, alpha=0.3)
ax5.set_xlim([0, 1])

# Plot hcov(3) relative error
ax6 = fig.add_subplot(gs[2, 1])
rel_error_3 = np.abs(hcov3_gvec - hcov3_vmec) / (np.abs(hcov3_vmec) + 1e-10)
ax6.semilogy(s_h, rel_error_3, 'g-', linewidth=2)
ax6.set_xlabel('s')
ax6.set_ylabel('Relative Error')
ax6.set_title('hcov(3) Relative Error')
ax6.grid(True, alpha=0.3)
ax6.set_xlim([0, 1])

plt.suptitle('hcov Components: VMEC vs GVEC Comparison\n(at θ=0, φ=0)', fontsize=14)
plt.tight_layout()
plt.savefig('hcov_1d_comparison.png', dpi=300, bbox_inches='tight')
plt.close()

# Create a summary plot showing all components together
fig2, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# Plot all components
ax1.plot(s_h, hcov1_vmec, 'b-', label='hcov(1) VMEC', linewidth=2)
ax1.plot(s_h, hcov1_gvec, 'b--', label='hcov(1) GVEC', linewidth=1.5, alpha=0.7)
ax1.plot(s_h, hcov2_vmec, 'r-', label='hcov(2) VMEC', linewidth=2)
ax1.plot(s_h, hcov2_gvec, 'r--', label='hcov(2) GVEC', linewidth=1.5, alpha=0.7)
ax1.plot(s_h, hcov3_vmec, 'g-', label='hcov(3) VMEC', linewidth=2)
ax1.plot(s_h, hcov3_gvec, 'g--', label='hcov(3) GVEC', linewidth=1.5, alpha=0.7)
ax1.set_ylabel('hcov components')
ax1.set_title('All hcov Components Comparison')
ax1.legend(ncol=2)
ax1.grid(True, alpha=0.3)
ax1.set_xlim([0, 1])

# Plot relative errors
rel_error_1 = np.abs(hcov1_gvec - hcov1_vmec) / (np.abs(hcov1_vmec) + 1e-10)
ax2.semilogy(s_h, rel_error_1, 'b-', label='hcov(1)', linewidth=2)
ax2.semilogy(s_h, rel_error_2, 'r-', label='hcov(2)', linewidth=2)
ax2.semilogy(s_h, rel_error_3, 'g-', label='hcov(3)', linewidth=2)
ax2.set_xlabel('s')
ax2.set_ylabel('Relative Error')
ax2.set_title('Relative Errors in hcov Components')
ax2.legend()
ax2.grid(True, alpha=0.3)
ax2.set_xlim([0, 1])
ax2.set_ylim([1e-6, 1e0])

plt.tight_layout()
plt.savefig('hcov_all_components.png', dpi=300, bbox_inches='tight')
plt.close()

print("hcov 1D comparison plots saved:")
print("  - hcov_1d_comparison.png")
print("  - hcov_all_components.png")