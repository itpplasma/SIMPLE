#!/usr/bin/env python3
"""
Plot results from existing test runs comparing VMEC internal Boozer vs BOOZXFORM
"""
import numpy as np
import matplotlib.pyplot as plt
import os

def load_confined_fraction(filename):
    """Load confined fraction data"""
    if os.path.exists(filename):
        data = np.loadtxt(filename)
        return data
    return None

def plot_comparison():
    """Plot comparison of existing test results"""
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('SIMPLE Field Type Comparison: Internal Boozer vs BOOZXFORM', fontsize=14)
    
    # Load data from test directories
    test_dirs = {
        'Internal Boozer': '/Users/ert/code/SIMPLE/test/booz_xform',
        'Test Run': '/Users/ert/code/SIMPLE/test_boozxform'
    }
    
    colors = {'Internal Boozer': 'blue', 'Test Run': 'red'}
    
    # Plot 1: Confined fraction vs time
    ax = axes[0, 0]
    for label, test_dir in test_dirs.items():
        conf_file = os.path.join(test_dir, 'confined_fraction.dat')
        if os.path.exists(conf_file):
            data = load_confined_fraction(conf_file)
            if data is not None and len(data) > 0:
                ax.plot(data[:, 0], data[:, 1], color=colors[label], 
                       label=label, linewidth=2)
    
    ax.set_xlabel('Time')
    ax.set_ylabel('Confined Fraction')
    ax.set_title('Confinement Evolution')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 2: Loss times
    ax = axes[0, 1]
    for label, test_dir in test_dirs.items():
        times_file = os.path.join(test_dir, 'times_lost.dat')
        if os.path.exists(times_file):
            try:
                times = np.loadtxt(times_file)
                if times.size > 0:
                    if times.ndim == 1:
                        # Single particle
                        ax.scatter([0], [times[0]], color=colors[label], 
                                 label=label, s=50)
                    else:
                        # Multiple particles
                        ax.scatter(range(len(times[:, 0])), times[:, 0], 
                                 color=colors[label], label=label, alpha=0.6, s=30)
            except:
                pass
    
    ax.set_xlabel('Particle Index')
    ax.set_ylabel('Loss Time')
    ax.set_title('Particle Loss Times')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 3: Initial conditions
    ax = axes[1, 0]
    for label, test_dir in test_dirs.items():
        start_file = os.path.join(test_dir, 'start.dat')
        if os.path.exists(start_file):
            try:
                start = np.loadtxt(start_file)
                if start.size > 0:
                    if start.ndim == 1:
                        # Single particle: s, R, Z, par_inv, perp_inv
                        ax.scatter([start[1]], [start[2]], color=colors[label], 
                                 label=f"{label} (s={start[0]:.2f})", s=100, marker='o')
                    else:
                        # Multiple particles
                        ax.scatter(start[:, 1], start[:, 2], color=colors[label], 
                                 label=label, alpha=0.6, s=30)
            except Exception as e:
                print(f"Error reading {start_file}: {e}")
    
    ax.set_xlabel('R [cm]')
    ax.set_ylabel('Z [cm]')
    ax.set_title('Initial Particle Positions')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.axis('equal')
    
    # Plot 4: Summary
    ax = axes[1, 1]
    ax.axis('off')
    
    summary_text = "Test Summary:\n\n"
    
    # Check what files exist
    for label, test_dir in test_dirs.items():
        summary_text += f"{label}:\n"
        
        # Check simple.in
        simple_in = os.path.join(test_dir, 'simple.in')
        if os.path.exists(simple_in):
            with open(simple_in, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if 'isw_field_type' in line:
                        summary_text += f"  Field type: {line.strip()}\n"
                    if 'boozxform_file' in line:
                        summary_text += f"  BOOZXFORM: {line.strip()}\n"
        
        # Check output files
        for fname in ['confined_fraction.dat', 'times_lost.dat', 'start.dat']:
            fpath = os.path.join(test_dir, fname)
            if os.path.exists(fpath):
                fsize = os.path.getsize(fpath)
                summary_text += f"  {fname}: {fsize} bytes\n"
        
        summary_text += "\n"
    
    ax.text(0.05, 0.95, summary_text, transform=ax.transAxes, 
            fontsize=9, verticalalignment='top', fontfamily='monospace')
    
    plt.tight_layout()
    return fig

def main():
    print("=== Plotting Existing SIMPLE Test Results ===")
    
    # Create comparison plots
    fig = plot_comparison()
    
    # Save figure
    output_file = 'field_type_comparison.png'
    fig.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\nSaved plot to {output_file}")
    
    # Also display
    plt.show()

if __name__ == '__main__':
    main()