#!/usr/bin/env python3
"""
Example: Plot energy loss vs J_perp from SIMPLE simulation output.

This script demonstrates how to visualize energy loss distribution
as a function of the perpendicular adiabatic invariant J_perp,
comparing runs with and without collisions.

The plot shows up to 4 curves:
1. COLL (actual): Energy from actual final velocity at loss (final_p^2)
2. NOCOLL (actual): Energy from actual final velocity (should be ~1)
3. COLL (theoretical): Energy from slowing-down curve lookup by loss time
4. NOCOLL (theoretical): Energy from slowing-down curve lookup

The theoretical curves use the bundled alpha particle slowing-down
curve (energyslow_aver.dat) to compute expected energy at loss time.

Usage:
    python plot_energy_loss.py /path/to/coll_run /path/to/nocoll_run
    python plot_energy_loss.py /path/to/coll_run  # COLL only

The output directories should contain:
    - times_lost.dat
    - confined_fraction.dat
    - start.dat (optional)

Requirements:
    - numpy
    - matplotlib
    - scipy
"""

import sys
from pathlib import Path

# Add parent directory to path for development
sys.path.insert(0, str(Path(__file__).parent.parent / "python"))

from pysimple.plotting import (
    load_loss_data,
    load_slowing_down_curve,
    plot_energy_loss_vs_jperp,
    plot_confined_fraction,
)


def main():
    if len(sys.argv) < 2:
        print("Usage: python plot_energy_loss.py <coll_dir> [nocoll_dir]")
        print("")
        print("Arguments:")
        print("  coll_dir    Directory with collision run output")
        print("  nocoll_dir  Directory with collisionless run output (optional)")
        sys.exit(1)

    coll_dir = Path(sys.argv[1])
    nocoll_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else None

    print(f"Loading COLL data from: {coll_dir}")
    data_coll = load_loss_data(coll_dir)

    data_nocoll = None
    if nocoll_dir:
        print(f"Loading NOCOLL data from: {nocoll_dir}")
        data_nocoll = load_loss_data(nocoll_dir)

    # Print summary
    print(f"\nCOLL run summary:")
    print(f"  Total particles: {data_coll.n_particles}")
    print(f"  Lost particles:  {data_coll.lost_mask.sum()}")
    print(f"  Confined:        {data_coll.confined_mask.sum()}")
    if data_coll.lost_mask.sum() > 0:
        fp = data_coll.final_p[data_coll.lost_mask]
        print(f"  final_p range:   {fp.min():.3f} - {fp.max():.3f}")

    if data_nocoll:
        print(f"\nNOCOLL run summary:")
        print(f"  Total particles: {data_nocoll.n_particles}")
        print(f"  Lost particles:  {data_nocoll.lost_mask.sum()}")
        print(f"  Confined:        {data_nocoll.confined_mask.sum()}")
        if data_nocoll.lost_mask.sum() > 0:
            fp = data_nocoll.final_p[data_nocoll.lost_mask]
            print(f"  final_p range:   {fp.min():.3f} - {fp.max():.3f}")

    # Load bundled slowing-down curve (no path needed)
    print("\nLoading bundled slowing-down curve...")
    sd_curve = load_slowing_down_curve()
    print(f"  Time range: {sd_curve[0][0]:.4f} - {sd_curve[0][-1]:.4f} s")
    print(f"  Energy range: {sd_curve[1][-1]:.4f} - {sd_curve[1][0]:.4f}")

    # Generate energy loss vs J_perp plot
    print("\nGenerating energy loss vs J_perp plot...")
    output_path = coll_dir / "energy_loss_vs_jperp.png"
    plot_energy_loss_vs_jperp(
        data_coll,
        data_nocoll,
        slowing_down_curve=sd_curve,
        title=coll_dir.name,
        output_path=output_path,
        show=False,
    )
    print(f"Saved: {output_path}")

    # Generate confined fraction plot with actual vs theoretical energy
    print("Generating confined fraction plot...")
    output_path = coll_dir / "confined_fraction_energy.png"
    plot_confined_fraction(
        data_coll,
        slowing_down_curve=sd_curve,
        title=f"{coll_dir.name} - Confined Fraction",
        output_path=output_path,
        show=False,
    )
    print(f"Saved: {output_path}")

    print("\nDone!")


if __name__ == "__main__":
    main()
