#!/usr/bin/env python3
"""
Example: Plot loss statistics from SIMPLE simulation output.

This script demonstrates how to use pysimple.plotting to visualize
particle loss statistics including:
- Confined fraction over time (particles AND energy)
- KDE density plot of loss time vs trapping parameter
- Energy loss distribution vs J_perp
- Starting position visualization

Usage:
    python plot_losses.py /path/to/simulation/output/
    python plot_losses.py  # uses current directory

The output directory should contain:
    - times_lost.dat
    - confined_fraction.dat
    - start.dat (optional but recommended)

Requirements:
    - numpy
    - matplotlib
    - scipy
"""

import sys
from pathlib import Path

# Add parent directory to path for development
sys.path.insert(0, str(Path(__file__).parent.parent / "python"))

# Import plotting module directly (does not require compiled Fortran backend)
import importlib.util
plotting_path = Path(__file__).parent.parent / "python" / "pysimple" / "plotting.py"
spec = importlib.util.spec_from_file_location("pysimple.plotting", plotting_path)
plotting = importlib.util.module_from_spec(spec)
sys.modules["pysimple.plotting"] = plotting
spec.loader.exec_module(plotting)

load_loss_data = plotting.load_loss_data
plot_loss_statistics = plotting.plot_loss_statistics
plot_confined_fraction = plotting.plot_confined_fraction
plot_kde_loss_density = plotting.plot_kde_loss_density
compute_energy_confined_fraction = plotting.compute_energy_confined_fraction


def main():
    # Get directory from command line or use current directory
    if len(sys.argv) > 1:
        data_dir = Path(sys.argv[1])
    else:
        data_dir = Path(".")

    print(f"Loading loss data from: {data_dir}")

    # Load all loss statistics data
    data = load_loss_data(data_dir)

    # Print summary
    print(f"\nData summary:")
    print(f"  Total particles: {data.n_particles}")
    print(f"  Lost particles:  {data.lost_mask.sum()}")
    print(f"  Confined:        {data.confined_mask.sum()}")
    print(f"  Skipped:         {data.skipped_mask.sum()}")
    print(f"  Trace time:      {data.trace_time:.3e} s")

    # Compute energy confined fraction
    t_energy, energy_frac = compute_energy_confined_fraction(data)
    final_particle_frac = data.confined_fraction[-1] if len(data.confined_fraction) > 0 else 0
    final_energy_frac = energy_frac[-1] if len(energy_frac) > 0 else 0

    print(f"\nFinal confinement:")
    print(f"  Particle confined fraction: {final_particle_frac:.3f}")
    print(f"  Energy confined fraction:   {final_energy_frac:.3f}")

    # Create comprehensive 4-panel figure
    print("\nGenerating loss statistics plot...")
    output_path = data_dir / "loss_statistics.png"
    fig = plot_loss_statistics(data, output_path=output_path, show=False)
    print(f"Saved: {output_path}")

    # Also create individual plots
    print("Generating confined fraction plot...")
    output_path = data_dir / "confined_fraction.png"
    plot_confined_fraction(data, output_path=output_path, show=False)
    print(f"Saved: {output_path}")

    print("Generating KDE density plot...")
    output_path = data_dir / "kde_loss_density.png"
    plot_kde_loss_density(data, output_path=output_path, show=False)
    print(f"Saved: {output_path}")

    print("\nDone!")


if __name__ == "__main__":
    main()
