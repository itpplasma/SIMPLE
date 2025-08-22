#!/usr/bin/env python3
"""
Simple script to compare orbits from VMEC vs BOOZXFORM field types
by running SIMPLE twice and plotting the results
"""
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os

def run_simple_with_field_type(field_type, label):
    """Run SIMPLE with specified field type and save orbit data"""
    
    # Create input file
    config = f"""&config
ntimstep = 1000
trace_time = 1d-4
sbeg = 0.5d0
ntestpart = 1
netcdffile = 'wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc'
boozxform_file = 'boozmn_LandremanPaul2021_QA_reactorScale_lowres_reference.nc'
isw_field_type = {field_type}
/"""
    
    with open('simple.in', 'w') as f:
        f.write(config)
    
    print(f"\nRunning SIMPLE with {label} (isw_field_type={field_type})...")
    
    # Run SIMPLE
    result = subprocess.run(['../build/simple.x'], 
                          capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error running SIMPLE with {label}:")
        print(result.stderr)
        return None
    
    # Read orbit data from fort.601 (or appropriate output file)
    orbit_file = 'fort.601'  # Poincare plot data
    if os.path.exists(orbit_file):
        try:
            data = np.loadtxt(orbit_file)
            # Save with unique name
            output_file = f'orbit_{label.lower().replace(" ", "_")}.dat'
            np.savetxt(output_file, data)
            print(f"  Saved orbit data to {output_file}")
            return data
        except:
            print(f"  Could not read orbit data from {orbit_file}")
            return None
    else:
        print(f"  No orbit output file {orbit_file} found")
        # Try to read from other possible output files
        for fortfile in ['fort.1001', 'fort.10001']:
            if os.path.exists(fortfile):
                try:
                    data = np.loadtxt(fortfile)
                    output_file = f'orbit_{label.lower().replace(" ", "_")}.dat'
                    np.savetxt(output_file, data)
                    print(f"  Found orbit data in {fortfile}, saved to {output_file}")
                    return data
                except:
                    pass
        return None

def plot_orbit_comparison(vmec_data, booz_data):
    """Create comparison plots of orbit data"""
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Orbit Comparison: VMEC vs BOOZXFORM Field Types', fontsize=14)
    
    # If we have Poincare data (R, Z at fixed toroidal angle)
    if vmec_data is not None and vmec_data.shape[1] >= 2:
        ax = axes[0, 0]
        if vmec_data.shape[1] >= 3:
            # Assume columns are: time/index, R, Z
            ax.plot(vmec_data[:, 1], vmec_data[:, 2], 'b.', 
                   markersize=3, label='VMEC', alpha=0.6)
        else:
            # Just two columns: R, Z
            ax.plot(vmec_data[:, 0], vmec_data[:, 1], 'b.', 
                   markersize=3, label='VMEC', alpha=0.6)
        ax.set_xlabel('R [cm]')
        ax.set_ylabel('Z [cm]')
        ax.set_title('Poincaré Section - VMEC')
        ax.axis('equal')
        ax.grid(True, alpha=0.3)
    
    if booz_data is not None and booz_data.shape[1] >= 2:
        ax = axes[0, 1]
        if booz_data.shape[1] >= 3:
            ax.plot(booz_data[:, 1], booz_data[:, 2], 'r.', 
                   markersize=3, label='BOOZXFORM', alpha=0.6)
        else:
            ax.plot(booz_data[:, 0], booz_data[:, 1], 'r.', 
                   markersize=3, label='BOOZXFORM', alpha=0.6)
        ax.set_xlabel('R [cm]')
        ax.set_ylabel('Z [cm]')
        ax.set_title('Poincaré Section - BOOZXFORM')
        ax.axis('equal')
        ax.grid(True, alpha=0.3)
    
    # Combined plot
    if vmec_data is not None and booz_data is not None:
        ax = axes[1, 0]
        if vmec_data.shape[1] >= 3 and booz_data.shape[1] >= 3:
            ax.plot(vmec_data[:, 1], vmec_data[:, 2], 'b.', 
                   markersize=3, label='VMEC', alpha=0.5)
            ax.plot(booz_data[:, 1], booz_data[:, 2], 'r.', 
                   markersize=3, label='BOOZXFORM', alpha=0.5)
        elif vmec_data.shape[1] >= 2 and booz_data.shape[1] >= 2:
            ax.plot(vmec_data[:, 0], vmec_data[:, 1], 'b.', 
                   markersize=3, label='VMEC', alpha=0.5)
            ax.plot(booz_data[:, 0], booz_data[:, 1], 'r.', 
                   markersize=3, label='BOOZXFORM', alpha=0.5)
        ax.set_xlabel('R [cm]')
        ax.set_ylabel('Z [cm]')
        ax.set_title('Overlay Comparison')
        ax.axis('equal')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    # Statistics
    ax = axes[1, 1]
    ax.axis('off')
    
    stats_text = "Orbit Statistics:\n\n"
    
    if vmec_data is not None:
        stats_text += f"VMEC field type:\n"
        stats_text += f"  Points: {len(vmec_data)}\n"
        if vmec_data.shape[1] >= 2:
            col_r = 1 if vmec_data.shape[1] >= 3 else 0
            col_z = 2 if vmec_data.shape[1] >= 3 else 1
            stats_text += f"  R range: [{vmec_data[:, col_r].min():.1f}, {vmec_data[:, col_r].max():.1f}]\n"
            stats_text += f"  Z range: [{vmec_data[:, col_z].min():.1f}, {vmec_data[:, col_z].max():.1f}]\n"
    else:
        stats_text += "VMEC: No data\n"
    
    stats_text += "\n"
    
    if booz_data is not None:
        stats_text += f"BOOZXFORM field type:\n"
        stats_text += f"  Points: {len(booz_data)}\n"
        if booz_data.shape[1] >= 2:
            col_r = 1 if booz_data.shape[1] >= 3 else 0
            col_z = 2 if booz_data.shape[1] >= 3 else 1
            stats_text += f"  R range: [{booz_data[:, col_r].min():.1f}, {booz_data[:, col_r].max():.1f}]\n"
            stats_text += f"  Z range: [{booz_data[:, col_z].min():.1f}, {booz_data[:, col_z].max():.1f}]\n"
    else:
        stats_text += "BOOZXFORM: No data\n"
    
    ax.text(0.1, 0.9, stats_text, transform=ax.transAxes, 
            fontsize=10, verticalalignment='top', fontfamily='monospace')
    
    plt.tight_layout()
    return fig

def main():
    print("=== SIMPLE Orbit Comparison: VMEC vs BOOZXFORM ===")
    
    # Check that files exist
    if not os.path.exists('wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc'):
        print("Error: VMEC file not found")
        return
    
    if not os.path.exists('boozmn_LandremanPaul2021_QA_reactorScale_lowres_reference.nc'):
        print("Error: BOOZXFORM file not found")
        return
    
    # Run with VMEC field type
    vmec_data = run_simple_with_field_type(1, "VMEC")
    
    # Run with BOOZXFORM field type
    booz_data = run_simple_with_field_type(5, "BOOZXFORM")
    
    # Create comparison plots
    if vmec_data is not None or booz_data is not None:
        print("\nCreating comparison plots...")
        fig = plot_orbit_comparison(vmec_data, booz_data)
        
        output_file = 'orbit_comparison.png'
        fig.savefig(output_file, dpi=150, bbox_inches='tight')
        print(f"Saved comparison plot to {output_file}")
        
        # Also show the plot
        plt.show()
    else:
        print("\nNo orbit data available for plotting")

if __name__ == '__main__':
    main()