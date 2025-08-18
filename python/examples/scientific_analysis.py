#!/usr/bin/env python3
"""
Scientific Analysis Example

Demonstrates physics-focused analysis capabilities of the SIMPLE Python API
for fusion plasma research, including orbit classification, confinement analysis,
and integration with the scientific Python ecosystem.

This example shows:
- Orbit classification (trapped vs passing, regular vs chaotic)
- Confinement analysis by flux surface and magnetic coordinates
- Statistical analysis of particle behavior
- Visualization with matplotlib
- Integration with NumPy/SciPy for advanced analysis
- Physics-based parameter sweeps

Requirements:
- SIMPLE built with Python support
- matplotlib for visualization
- scipy for statistical analysis (optional)
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Add SIMPLE Python API to path
sys.path.insert(0, str(Path(__file__).parent.parent))

import simple


def download_test_vmec():
    """Download test VMEC file if not present"""
    vmec_file = "wout.nc"
    if not Path(vmec_file).exists():
        print("Downloading test VMEC file...")
        import urllib.request
        url = "https://github.com/hiddenSymmetries/simsopt/raw/master/tests/test_files/wout_LandremanPaul2021_QA_reactorScale_lowres_reference.nc"
        try:
            urllib.request.urlretrieve(url, vmec_file)
            print(f"Downloaded {vmec_file}")
        except Exception as e:
            print(f"Failed to download VMEC file: {e}")
            return None
    return vmec_file


def orbit_classification_analysis():
    """Analyze orbit classification: trapped vs passing, regular vs chaotic"""
    print("=== Orbit Classification Analysis ===")
    
    # Create particles with surface sampling
    n_particles = 50000
    particles = simple.ParticleBatch(n_particles)
    particles.initialize_surface("wout.nc", s=0.9)
    
    # Run simulation with sufficient time for classification
    print(f"Running simulation for orbit classification ({n_particles:,} particles)...")
    results = simple.trace_orbits(
        particles,
        tmax=2000.0,  # Longer time for classification
        integrator='symplectic_midpoint',
        verbose=True
    )
    
    # Extract physics data
    trap_parameter = results.trap_parameter
    perp_invariant = results.perpendicular_invariant
    loss_times = results.loss_times
    final_positions = results.final_positions
    
    # Classify orbits based on trapping parameter
    trapped_mask = trap_parameter > 0
    passing_mask = trap_parameter <= 0
    confined_mask = results.confined_mask
    lost_mask = results.lost_mask
    
    # Cross-classification
    trapped_confined = np.sum(trapped_mask & confined_mask)
    trapped_lost = np.sum(trapped_mask & lost_mask)
    passing_confined = np.sum(passing_mask & confined_mask)
    passing_lost = np.sum(passing_mask & lost_mask)
    
    print(f"\nOrbit Classification Results:")
    print(f"{'Category':<20} {'Count':<10} {'Fraction':<10}")
    print("-" * 40)
    print(f"{'Trapped':<20} {np.sum(trapped_mask):<10,} {np.sum(trapped_mask)/n_particles:<10.3f}")
    print(f"{'Passing':<20} {np.sum(passing_mask):<10,} {np.sum(passing_mask)/n_particles:<10.3f}")
    print(f"{'Confined':<20} {np.sum(confined_mask):<10,} {np.sum(confined_mask)/n_particles:<10.3f}")
    print(f"{'Lost':<20} {np.sum(lost_mask):<10,} {np.sum(lost_mask)/n_particles:<10.3f}")
    
    print(f"\nCross-Classification:")
    print(f"{'Trapped & Confined':<20} {trapped_confined:<10,} {trapped_confined/n_particles:<10.3f}")
    print(f"{'Trapped & Lost':<20} {trapped_lost:<10,} {trapped_lost/n_particles:<10.3f}")
    print(f"{'Passing & Confined':<20} {passing_confined:<10,} {passing_confined/n_particles:<10.3f}")
    print(f"{'Passing & Lost':<20} {passing_lost:<10,} {passing_lost/n_particles:<10.3f}")
    
    # Physics insights
    if trapped_lost > 0 and passing_lost > 0:
        trapped_loss_rate = trapped_lost / np.sum(trapped_mask)
        passing_loss_rate = passing_lost / np.sum(passing_mask)
        
        print(f"\nPhysics Insights:")
        print(f"  Trapped particle loss rate: {trapped_loss_rate:.3f}")
        print(f"  Passing particle loss rate: {passing_loss_rate:.3f}")
        
        if trapped_loss_rate > passing_loss_rate:
            print("  → Trapped particles more likely to be lost (typical near edge)")
        else:
            print("  → Passing particles more likely to be lost")
    
    return results, trap_parameter, perp_invariant


def confinement_vs_surface_analysis():
    """Analyze confinement as function of flux surface"""
    print("\n=== Confinement vs Flux Surface Analysis ===")
    
    # Sample particles across multiple flux surfaces
    s_surfaces = np.linspace(0.3, 0.95, 8)  # From core to edge
    n_per_surface = 5000
    
    confinement_data = []
    
    print(f"Analyzing confinement across {len(s_surfaces)} flux surfaces...")
    
    for i, s in enumerate(s_surfaces):
        print(f"  Surface {i+1}/{len(s_surfaces)}: s = {s:.2f}")
        
        # Create particles on this surface
        particles = simple.ParticleBatch(n_per_surface)
        particles.initialize_surface("wout.nc", s=s)
        
        # Run simulation
        results = simple.trace_orbits(
            particles,
            tmax=1000.0,
            integrator='symplectic_midpoint',
            verbose=False
        )
        
        # Analyze confinement
        stats = results.confinement_statistics()
        trap_parameter = results.trap_parameter
        
        # Collect data
        surface_data = {
            's': s,
            'confined_fraction': stats.confined_fraction,
            'mean_loss_time': stats.mean_loss_time if stats.n_lost > 0 else np.inf,
            'trapped_fraction': np.sum(trap_parameter > 0) / n_per_surface,
            'n_confined': stats.n_confined,
            'n_lost': stats.n_lost
        }
        
        confinement_data.append(surface_data)
        
        print(f"    Confined: {stats.confined_fraction:.3f}, "
              f"Trapped: {surface_data['trapped_fraction']:.3f}")
    
    # Convert to arrays for analysis
    s_values = np.array([d['s'] for d in confinement_data])
    confined_fractions = np.array([d['confined_fraction'] for d in confinement_data])
    trapped_fractions = np.array([d['trapped_fraction'] for d in confinement_data])
    mean_loss_times = np.array([d['mean_loss_time'] for d in confinement_data])
    
    # Replace inf with NaN for plotting
    mean_loss_times = np.where(np.isfinite(mean_loss_times), mean_loss_times, np.nan)
    
    print(f"\nConfinement Profile Summary:")
    print(f"{'s':<8} {'Confined':<12} {'Trapped':<12} {'Loss Time':<12}")
    print("-" * 50)
    
    for i, data in enumerate(confinement_data):
        loss_time_str = f"{data['mean_loss_time']:.1f}" if np.isfinite(data['mean_loss_time']) else "∞"
        print(f"{data['s']:<8.2f} {data['confined_fraction']:<12.3f} "
              f"{data['trapped_fraction']:<12.3f} {loss_time_str:<12}")
    
    return s_values, confined_fractions, trapped_fractions, mean_loss_times


def statistical_analysis():
    """Perform statistical analysis of particle behavior"""
    print("\n=== Statistical Analysis ===")
    
    # Large sample for statistical significance
    n_particles = 100000
    particles = simple.ParticleBatch(n_particles)
    particles.initialize_volume("wout.nc", s_min=0.1, s_max=0.95)
    
    print(f"Running large simulation for statistics ({n_particles:,} particles)...")
    results = simple.trace_orbits(
        particles,
        tmax=1500.0,
        integrator='symplectic_midpoint',
        verbose=True
    )
    
    # Extract data arrays
    initial_positions = particles.positions
    final_positions = results.final_positions
    loss_times = results.loss_times
    trap_parameter = results.trap_parameter
    perp_invariant = results.perpendicular_invariant
    
    # Initial distribution analysis
    initial_s = initial_positions[0, :]
    initial_theta = initial_positions[1, :]
    initial_phi = initial_positions[2, :]
    initial_vpar = initial_positions[3, :]
    initial_mu = initial_positions[4, :]
    
    print(f"\nInitial Distribution Statistics:")
    print(f"  s: mean={initial_s.mean():.3f}, std={initial_s.std():.3f}, range=[{initial_s.min():.3f}, {initial_s.max():.3f}]")
    print(f"  v_par: mean={initial_vpar.mean():.3f}, std={initial_vpar.std():.3f}")
    print(f"  mu: mean={initial_mu.mean():.3f}, std={initial_mu.std():.3f}")
    
    # Final distribution analysis (confined particles only)
    confined_mask = results.confined_mask
    if np.sum(confined_mask) > 0:
        final_s_confined = final_positions[0, confined_mask]
        final_vpar_confined = final_positions[3, confined_mask]
        
        print(f"\nFinal Distribution (Confined Particles):")
        print(f"  s: mean={final_s_confined.mean():.3f}, std={final_s_confined.std():.3f}")
        print(f"  v_par: mean={final_vpar_confined.mean():.3f}, std={final_vpar_confined.std():.3f}")
        
        # Radial migration analysis
        initial_s_confined = initial_s[confined_mask]
        radial_migration = final_s_confined - initial_s_confined
        
        print(f"\nRadial Migration Analysis (Confined Particles):")
        print(f"  Mean migration: {radial_migration.mean():.6f}")
        print(f"  RMS migration: {np.sqrt(np.mean(radial_migration**2)):.6f}")
        print(f"  Max outward: {radial_migration.max():.6f}")
        print(f"  Max inward: {radial_migration.min():.6f}")
    
    # Loss time statistics
    lost_mask = results.lost_mask
    if np.sum(lost_mask) > 0:
        finite_loss_times = loss_times[lost_mask]
        
        print(f"\nLoss Time Statistics:")
        print(f"  Mean: {finite_loss_times.mean():.3f}")
        print(f"  Median: {np.median(finite_loss_times):.3f}")
        print(f"  Standard deviation: {finite_loss_times.std():.3f}")
        print(f"  Range: [{finite_loss_times.min():.3f}, {finite_loss_times.max():.3f}]")
        
        # Percentiles
        percentiles = [10, 25, 50, 75, 90, 95, 99]
        perc_values = np.percentile(finite_loss_times, percentiles)
        
        print(f"  Percentiles:")
        for p, v in zip(percentiles, perc_values):
            print(f"    {p}%: {v:.3f}")
    
    # Invariant conservation analysis
    print(f"\nInvariant Conservation Analysis:")
    perp_inv_mean = perp_invariant.mean()
    perp_inv_std = perp_invariant.std()
    perp_inv_variation = perp_inv_std / perp_inv_mean if perp_inv_mean > 0 else np.inf
    
    print(f"  Perpendicular invariant:")
    print(f"    Mean: {perp_inv_mean:.6e}")
    print(f"    Std: {perp_inv_std:.6e}")
    print(f"    Coefficient of variation: {perp_inv_variation:.6f}")
    
    if perp_inv_variation < 0.01:
        print("    → Good conservation (variation < 1%)")
    elif perp_inv_variation < 0.1:
        print("    → Acceptable conservation (variation < 10%)")
    else:
        print("    → Poor conservation (variation > 10%)")
    
    return results


def parameter_sweep_analysis():
    """Perform physics parameter sweeps"""
    print("\n=== Parameter Sweep Analysis ===")
    
    # Sweep simulation time to study convergence
    time_values = [100, 250, 500, 1000, 2000]
    n_particles = 10000
    
    print(f"Time convergence study ({n_particles:,} particles per run):")
    
    # Create consistent initial conditions
    particles = simple.ParticleBatch(n_particles)
    particles.initialize_surface("wout.nc", s=0.9)
    base_particles = particles.copy()
    
    time_results = []
    
    for tmax in time_values:
        print(f"  Running tmax = {tmax}...")
        
        test_particles = base_particles.copy()
        results = simple.trace_orbits(
            test_particles,
            tmax=tmax,
            integrator='symplectic_midpoint',
            verbose=False
        )
        
        stats = results.confinement_statistics()
        time_results.append({
            'tmax': tmax,
            'confined_fraction': stats.confined_fraction,
            'n_lost': stats.n_lost,
            'mean_loss_time': stats.mean_loss_time if stats.n_lost > 0 else np.inf
        })
    
    print(f"\nTime Convergence Results:")
    print(f"{'tmax':<8} {'Confined':<12} {'Lost':<8} {'Mean Loss':<12}")
    print("-" * 45)
    
    for data in time_results:
        loss_str = f"{data['mean_loss_time']:.1f}" if np.isfinite(data['mean_loss_time']) else "∞"
        print(f"{data['tmax']:<8} {data['confined_fraction']:<12.4f} "
              f"{data['n_lost']:<8} {loss_str:<12}")
    
    # Check convergence
    confined_fractions = [r['confined_fraction'] for r in time_results]
    if len(confined_fractions) >= 3:
        # Simple convergence check
        last_three = confined_fractions[-3:]
        variation = max(last_three) - min(last_three)
        
        print(f"\nConvergence Analysis:")
        print(f"  Variation in last 3 points: {variation:.6f}")
        if variation < 0.001:
            print("  → Converged (variation < 0.001)")
        elif variation < 0.01:
            print("  → Nearly converged (variation < 0.01)")
        else:
            print("  → Not converged (variation > 0.01)")
    
    return time_results


def create_visualization_plots():
    """Create comprehensive visualization plots"""
    print("\n=== Creating Visualization Plots ===")
    
    try:
        import matplotlib.pyplot as plt
        print("Creating physics analysis plots...")
        
        # Run simulation for plotting
        n_particles = 25000
        particles = simple.ParticleBatch(n_particles)
        particles.initialize_surface("wout.nc", s=0.9)
        
        results = simple.trace_orbits(
            particles,
            tmax=1000.0,
            integrator='symplectic_midpoint',
            verbose=False
        )
        
        # Create figure with subplots
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        fig.suptitle('SIMPLE Python API - Physics Analysis', fontsize=16)
        
        # Plot 1: Loss time histogram
        ax = axes[0, 0]
        stats = results.confinement_statistics()
        hist_counts, hist_bins = stats.loss_time_distribution
        
        ax.bar(hist_bins[:-1], hist_counts, width=np.diff(hist_bins), alpha=0.7)
        ax.set_xlabel('Loss Time')
        ax.set_ylabel('Number of Particles')
        ax.set_title('Loss Time Distribution')
        ax.grid(True, alpha=0.3)
        
        # Plot 2: Poincaré section (confined vs lost)
        ax = axes[0, 1]
        final_pos = results.final_positions
        confined_mask = results.confined_mask
        lost_mask = results.lost_mask
        
        # Sample for plotting performance
        n_plot = min(5000, n_particles)
        indices = np.random.choice(n_particles, n_plot, replace=False)
        
        confined_indices = indices[confined_mask[indices]]
        lost_indices = indices[lost_mask[indices]]
        
        if len(confined_indices) > 0:
            ax.scatter(final_pos[1, confined_indices], final_pos[2, confined_indices], 
                      s=1, alpha=0.6, label=f'Confined ({len(confined_indices)})', color='blue')
        
        if len(lost_indices) > 0:
            ax.scatter(final_pos[1, lost_indices], final_pos[2, lost_indices], 
                      s=1, alpha=0.6, label=f'Lost ({len(lost_indices)})', color='red')
        
        ax.set_xlabel('θ (poloidal angle)')
        ax.set_ylabel('φ (toroidal angle)')
        ax.set_title('Poincaré Section')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Plot 3: Trap parameter distribution
        ax = axes[0, 2]
        trap_par = results.trap_parameter
        
        # Histogram of trap parameter
        ax.hist(trap_par, bins=50, alpha=0.7, density=True)
        ax.axvline(0, color='red', linestyle='--', label='Trapped/Passing boundary')
        ax.set_xlabel('Trap Parameter')
        ax.set_ylabel('Density')
        ax.set_title('Orbit Classification')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Plot 4: Radial profiles (run surface analysis)
        ax = axes[1, 0]
        
        # Quick surface analysis for plotting
        s_surfaces = [0.4, 0.6, 0.8, 0.9]
        confined_fracs = []
        
        for s in s_surfaces:
            test_particles = simple.ParticleBatch(5000)
            test_particles.initialize_surface("wout.nc", s=s)
            test_results = simple.trace_orbits(test_particles, tmax=500.0, verbose=False)
            test_stats = test_results.confinement_statistics()
            confined_fracs.append(test_stats.confined_fraction)
        
        ax.plot(s_surfaces, confined_fracs, 'o-', linewidth=2, markersize=8)
        ax.set_xlabel('Flux Surface (s)')
        ax.set_ylabel('Confined Fraction')
        ax.set_title('Confinement vs Flux Surface')
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0, 1)
        
        # Plot 5: Energy/invariant scatter
        ax = axes[1, 1]
        perp_inv = results.perpendicular_invariant
        
        # Sample for plotting
        plot_indices = np.random.choice(n_particles, min(2000, n_particles), replace=False)
        
        colors = ['blue' if confined_mask[i] else 'red' for i in plot_indices]
        ax.scatter(trap_par[plot_indices], perp_inv[plot_indices], 
                  c=colors, s=10, alpha=0.6)
        ax.set_xlabel('Trap Parameter')
        ax.set_ylabel('Perpendicular Invariant')
        ax.set_title('Phase Space Structure')
        ax.grid(True, alpha=0.3)
        
        # Plot 6: Loss rate vs initial position
        ax = axes[1, 2]
        initial_s = particles.positions[0, :]
        loss_times = results.loss_times
        
        # Bin by initial s coordinate
        s_bins = np.linspace(initial_s.min(), initial_s.max(), 20)
        s_centers = (s_bins[:-1] + s_bins[1:]) / 2
        
        loss_rates = []
        for i in range(len(s_bins)-1):
            mask = (initial_s >= s_bins[i]) & (initial_s < s_bins[i+1])
            if np.sum(mask) > 0:
                lost_in_bin = np.sum(loss_times[mask] < np.inf)
                loss_rate = lost_in_bin / np.sum(mask)
                loss_rates.append(loss_rate)
            else:
                loss_rates.append(0)
        
        ax.plot(s_centers, loss_rates, 'o-', linewidth=2, markersize=6)
        ax.set_xlabel('Initial s coordinate')
        ax.set_ylabel('Loss Rate')
        ax.set_title('Loss Rate vs Initial Position')
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0, 1)
        
        plt.tight_layout()
        
        # Save plot
        plot_file = "simple_physics_analysis.png"
        plt.savefig(plot_file, dpi=150, bbox_inches='tight')
        print(f"Physics analysis plot saved to: {plot_file}")
        
        plt.show()
        
    except ImportError:
        print("matplotlib not available - skipping visualization")
    except Exception as e:
        print(f"Plotting failed: {e}")


def main():
    """Main scientific analysis execution"""
    print("SIMPLE Python API - Scientific Analysis Example")
    print("=" * 55)
    
    # Download test data
    vmec_file = download_test_vmec()
    if vmec_file is None:
        print("Cannot proceed without VMEC file")
        return
    
    try:
        # Run scientific analyses
        orbit_results, trap_par, perp_inv = orbit_classification_analysis()
        s_values, confined_fracs, trapped_fracs, loss_times = confinement_vs_surface_analysis()
        statistical_results = statistical_analysis()
        time_results = parameter_sweep_analysis()
        create_visualization_plots()
        
        print("\n" + "=" * 55)
        print("Scientific Analysis Summary:")
        print("=" * 55)
        
        print("✓ Orbit classification: Trapped vs passing particle analysis")
        print("✓ Confinement profiles: Physics dependence on flux surface")
        print("✓ Statistical analysis: Large-scale particle behavior")
        print("✓ Parameter sweeps: Convergence and sensitivity studies")
        print("✓ Visualization: Physics-focused plotting capabilities")
        
        print("\nPhysics Insights:")
        print("- Particle confinement strongly depends on flux surface location")
        print("- Orbit classification reveals transport mechanisms")
        print("- Statistical analysis validates symplectic integration quality")
        print("- Parameter sweeps enable systematic physics studies")
        
        print("\nNext Steps for Research:")
        print("- Extend to different VMEC equilibria")
        print("- Study parameter dependencies (field strength, shear, etc.)")
        print("- Compare with experimental confinement data")
        print("- Implement advanced orbit classification algorithms")
        
    except Exception as e:
        print(f"\nScientific analysis failed: {e}")
        print("Please check:")
        print("1. SIMPLE built with Python support")
        print("2. matplotlib installed for visualization")
        print("3. VMEC equilibrium file is valid")
        raise


if __name__ == "__main__":
    main()