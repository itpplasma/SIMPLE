# Wall Heat Load Physics

This document describes the physics behind the wall heat load calculation in `pysimple.engineering`.

## D-T Fusion Energy Balance

Each deuterium-tritium fusion reaction produces 17.6 MeV:

$$\text{D} + \text{T} \rightarrow \text{He}^{4} (3.5\,\text{MeV}) + n (14.1\,\text{MeV})$$

- **Alpha particle** ($\alpha$): 3.5 MeV kinetic energy
- **Neutron**: 14.1 MeV (escapes plasma, heats blanket)

The alpha power is therefore:

$$P_\alpha = \frac{P_\text{fusion}}{5}$$

For a 3 GW fusion reactor: $P_\alpha \approx 600$ MW.

## Alpha Particle Slowing Down

Alphas are born at $E_0 = 3.5$ MeV and slow down via Coulomb collisions with background electrons and ions. The characteristic slowing-down time is:

$$\tau_s = \frac{3(2\pi)^{3/2} \epsilon_0^2 m_\alpha}{e^4 Z_\alpha^2 n_e \ln\Lambda} \left(\frac{T_e}{m_e}\right)^{3/2}$$

For typical reactor parameters ($n_e = 10^{20}$ m$^{-3}$, $T_e = 10$ keV), $\tau_s \approx 0.1-1$ s.

The implemented operator includes momentum drag, momentum diffusion and
pitch scattering. Its energy histories depend on the supplied temperature
and density profiles; the characteristic time is a guide for choosing the
trace duration, not an exact stopping law for each marker.

## Collisions with an STL wall

Build with `SIMPLE_ENABLE_CGAL=ON`. Set `swcoll = .true.` and supply both
`wall_input` (an STL file) and `coord_input` (a chartmap covering the plasma
and the wall). The field must remain valid throughout that chart, for example
by splining a coil field onto an extended chartmap. Set `wall_units = 'm'` or
`'mm'` to match the STL. A wall outside the LCFS needs an extended chart;
reaching the chart boundary alone does not establish an STL hit.

The same input file contains the background plasma profiles, in eV and m$^{-3}$:

```fortran
&config
  swcoll = .true.
  field_input = 'coils.simple'
  coord_input = 'wall.chartmap.nc'
  wall_input = 'wall.stl'
  wall_units = 'm'
  isw_field_type = 5  ! Reference-coordinate field
  integmode = 0      ! RK guiding-center integration
  facE_al = 1        ! Birth energy 3.5 MeV
  contr_pp = -1d10   ! Trace strongly passing markers too
  trace_time = 0.1
  deterministic = .true.
  ran_seed = 541
  output_results_netcdf = .true.
/
&profiles
  profile_type = 'power_series'
  Te_scale = 1d4
  Ti1_scale = 1d4
  Ti2_scale = 1d4
  ni1_scale = 5d19
  ni2_scale = 5d19
/
```

The default constant coefficients make these flat profiles. Choose
`trace_time` relative to the slowing-down time of the supplied profiles:
`trace_time << tau_s` measures prompt losses with little slowing down;
times comparable to `tau_s` resolve losses at reduced energy. Check the
birth momentum slowing rate (`params::slowrate`, in s$^{-1}$, computed at
collision initialization); its inverse
is a local momentum relaxation scale, not a universal stopping time. Energy
diffusion and the spatial profiles preclude a single deterministic energy law
for all markers. Resolve collision kicks with the microstep size and check
convergence when changing the density, temperature, or trace time.

Each microstep first advances the guiding center and checks the Cartesian
segment for its first STL intersection. A hit retains momentum and pitch
from that drift, including collision kicks from preceding microsteps. The
code does not apply the full next collision kick after an already detected
hit, or a fractional kick on the final segment. Thus hit energy has the
time resolution of the drift/kick splitting. `times_lost` is interpolated
along the intersecting segment. `wall_hit_angle_rad` is the angle between
the segment and the wall normal (0 is normal incidence), not a full-orbit
velocity incidence angle.

Use `wall_hit == 1` to select actual STL intersections; `class_lost` is a
classifier flag and does not select STL hits in ordinary tracing. Heat-load
weights use the recorded
`zend(4)**2 = E_hit/E0`, with `E0 = 3.5 MeV / facE_al`. Energy transferred
to the background before impact is not deposited in the wall bin. Thermal
diffusion can increase individual energies, so `p < 1` is an ensemble
slowing-down expectation rather than a strict bound for every marker.

The current `WallHeatMap.from_netcdf` reader selects `class_lost`, so ordinary
STL tracing results need direct selection with `wall_hit` for heat-load
post-processing. The regression CSV supplies the hit energy weights directly.

The registered `test_collisional_stl_wall` case uses 256 3.5 MeV alphas,
the NCSX coil fixture, and an STL triangulated from the boundary of
the map2disc chart extended by 0.1 m, contracted 5% toward the chart axis
to leave a field-domain buffer around the planar facets. It uses enhanced flat
densities and 1 keV temperatures to resolve slowing down in a short trace;
this is a regression fixture, not a reactor prediction. With Python tools,
map2disc, NumPy and netCDF4 available, run (the exercised map2disc environment
used NumPy 2.2.6 and SciPy 1.16.3):

```sh
CMAKE_BUILD_PARALLEL_LEVEL=8 make test CONFIG=Fast \
  FLAGS='-DSIMPLE_ENABLE_CGAL=ON -DENABLE_PYTHON_INTERFACE=OFF -DSIMPLE_ENABLE_PYTHON_TOOLS=ON' \
  TEST=test_collisional_stl_wall
```

CTest generates the chartmap fixture automatically. The test checks the
collisionless momentum invariant, the negligible-density limit, reduced
energy on delayed wall hits, and identical particle results in two seeded
one-thread runs. It writes `build/test/tests/collisional_stl_wall/wall_hits.csv`
with mode, marker, hit energy fraction, hit time in seconds, and incidence
angle in radians, alongside per-run `results.nc` and `run.log`. The fixed
thread count is `OMP_NUM_THREADS=1`; multithread stream independence and
scheduling-dependent reproducibility remain tracked in #42.

## Monte Carlo Heat Flux Calculation

SIMPLE traces $N$ test particles representing the alpha population. Each particle carries a fraction of total alpha power:

$$\text{Power per MC particle} = \frac{P_\alpha}{N}$$

### Energy Weighting

Particles that escape early (before significant slowing down) carry more energy than those that escape late. The normalized velocity at wall impact is:

$$p = \frac{v}{v_0} = \sqrt{\frac{E}{E_0}}$$

where $v_0$ is the birth velocity. The energy carried by particle $i$ at wall impact is:

$$E_i = p_i^2 \cdot E_0$$

Without collisions: $p = 1$ for all particles.
With collisions: $p < 1$ for particles that slowed down before escaping.

### Heat Flux Calculation

For a wall bin with area $A_\text{bin}$, the energy-weighted heat flux is:

$$q_\text{bin} = \frac{1}{A_\text{bin}} \sum_{i \in \text{bin}} p_i^2 \cdot \frac{P_\alpha}{N}$$

Or equivalently:

$$q_\text{bin} = \frac{\sum_{i \in \text{bin}} p_i^2}{N} \cdot \frac{P_\alpha}{A_\text{bin}}$$

### Energy Loss Fraction

The **particle loss fraction** is:

$$f_\text{particle} = \frac{n_\text{lost}}{N}$$

The **energy loss fraction** accounts for slowing down:

$$f_\text{energy} = \frac{1}{N} \sum_{i \in \text{lost}} p_i^2$$

For an ensemble dominated by slowing down, typically
$f_\text{energy} \leq f_\text{particle}$. Stochastic energy diffusion does not
impose $p_i^2 \leq 1$ on every individual marker.

The lost power is:

$$P_\text{lost} = f_\text{energy} \cdot P_\alpha$$

## Typical Values

| Configuration | Peak Heat Flux | Reference |
|--------------|----------------|-----------|
| Infinity Two | ~2.5 MW/m$^2$ | Albert et al. (2024) |
| ARIES-CS | 5-18 MW/m$^2$ | Ku & Boozer (2008) |
| W7-X (experiments) | 0.1-1 MW/m$^2$ | Lazerson et al. (2021) |
| Material limits | 5-10 MW/m$^2$ | Engineering constraint |

## Wall Area Calculation

The heat flux requires accurate wall surface area. Two methods are available:

### Simple Torus Approximation

For quick estimates, the code defaults to a simple torus:

$$A = 4\pi^2 R_0 a$$

where $R_0$ is major radius and $a$ is minor radius. This underestimates stellarator wall area by 1.5-2x due to 3D shaping.

### Actual Geometry from Chartmap

For accurate results, provide a chartmap file with the actual wall geometry:

$$A = n_{fp} \iint \left|\frac{\partial \mathbf{r}}{\partial \theta} \times \frac{\partial \mathbf{r}}{\partial \zeta}\right| d\theta\, d\zeta$$

The chartmap contains Cartesian coordinates $(x, y, z)$ on a $(\rho, \theta, \zeta)$ grid. The wall surface is at $\rho = 1$.

### Resolution Matching

When using a chartmap, the heat map bin resolution should not exceed the chartmap grid resolution. If the heat map has more bins than chartmap grid points, many bins will have zero computed area, leading to inaccurate flux calculations.

The code issues a warning if the sum of bin areas is less than 90% of the total wall area:

```
UserWarning: Bin areas sum (X m^2) is less than 90% of wall area (Y m^2).
Heat map resolution may exceed chartmap resolution, causing empty bins.
Consider reducing n_theta/n_zeta or using a finer chartmap.
```

**Recommendation**: Use `n_theta` ≤ chartmap ntheta and `n_zeta` ≤ chartmap nzeta × nfp.

## Guiding Center Limitation

**Important**: SIMPLE tracks guiding center orbits, not full particle orbits. The reported wall positions are guiding center positions, not actual wall impact points.

For 3.5 MeV alpha particles in a ~5 T field:

$$\rho_\alpha = \frac{m_\alpha v_\perp}{Z_\alpha e B} \approx 5-10\,\text{cm}$$

This means:
- Heat load patterns are smeared by ~10 cm
- Fine-scale wall features are not resolved
- Peak flux may be underestimated due to averaging

For high-resolution wall load studies (e.g., limiter tile design), finite Larmor radius effects should be added by:
1. Adding a gyroradius offset in the direction perpendicular to B
2. Or using a full-orbit code instead of guiding center

## Code Implementation

The `WallHeatMap` class in `pysimple.engineering` implements this calculation:

```python
from pysimple.engineering import WallHeatMap, compute_wall_area_from_chartmap

# With actual wall geometry from chartmap (recommended)
heat_map = WallHeatMap.from_netcdf(
    "results.nc",
    total_alpha_power_MW=600.0,  # 3 GW fusion -> 600 MW alpha
    chartmap_file="wout.chartmap.nc",  # Actual 3D wall geometry
)

# Or with simple torus approximation (for quick estimates)
heat_map_simple = WallHeatMap.from_netcdf(
    "results.nc",
    total_alpha_power_MW=600.0,
    major_radius=1000.0,  # cm
    minor_radius=100.0,   # cm
)

print(f"Particle loss fraction: {heat_map.loss_fraction:.1%}")
print(f"Energy loss fraction: {heat_map.energy_loss_fraction:.1%}")
print(f"Lost power: {heat_map.lost_power:.1f} MW")
print(f"Peak flux: {heat_map.peak_flux:.2f} MW/m^2")
print(f"Wall area: {heat_map.wall_area:.1f} m^2")
```

## Data Format

The NetCDF results file from SIMPLE contains:

| Variable | Shape | Description |
|----------|-------|-------------|
| `zend` | (5, N) | Final reference phase space: radial coordinate (rho for chartmaps), theta, zeta, p=v/v0, lambda |
| `xend_cart` | (3, N) | Final Cartesian position (x, y, z) |
| `class_lost` | (N,) | Classifier loss flag; use wall_hit to select STL intersections |
| `wall_hit` | (N,) | 1 for a confirmed STL intersection |
| `wall_hit_angle_rad` | (N,) | Segment incidence angle relative to wall normal (rad) |
| `times_lost` | (N,) | Time of wall impact (s) |

The velocity ratio `zend[3]` = $p = v/v_0$ is crucial for energy weighting.

## References

1. **Lazerson, S.A., et al.** (2021). "Modeling of energetic particle transport in optimized stellarators." *Plasma Physics and Controlled Fusion* 63, 125033.
   https://doi.org/10.1088/1361-6587/ac2b38

2. **Albert, C.G., et al.** (2024). "Alpha particle confinement in Infinity Two." *Journal of Plasma Physics*.
   https://doi.org/10.1017/S0022377824000163

3. **Ku, L.P. & Boozer, A.H.** (2008). "Stellarator alpha particle loss calculations for ARIES-CS." *Fusion Science and Technology* 54, 673-693.
   https://doi.org/10.13182/FST08-A1891

4. **Drevlak, M., et al.** (2014). "Fast particle confinement with optimized coil currents in the W7-X stellarator." *Nuclear Fusion* 54, 073002.
   https://doi.org/10.1088/0029-5515/54/7/073002

5. **BEAMS3D code documentation**: https://princetonuniversity.github.io/STELLOPT/BEAMS3D
