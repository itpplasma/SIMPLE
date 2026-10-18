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

The energy evolution follows:

$$E(t) = E_0 \left(1 - t/\tau_s\right)^{2/3} \quad \text{for } t < \tau_s$$

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

Note: $f_\text{energy} \leq f_\text{particle}$ because $p_i^2 \leq 1$.

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
| `zend` | (5, N) | Final phase space: s, theta, zeta, p=v/v0, lambda |
| `xend_cart` | (3, N) | Final Cartesian position (x, y, z) |
| `class_lost` | (N,) | 1 if particle hit wall, 0 otherwise |
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

## Heat loads on STL triangles

`pysimple.wall_heat_load.TriangleHeatLoad` uses the engineering mesh directly.
Pure post-processing imports work without the compiled SIMPLE backend; calling
the tracing API still requires it.
`pysimple.engineering.WallHeatMap` remains available for flux-coordinate bins.
Build SIMPLE with `SIMPLE_ENABLE_CGAL=ON`, set `wall_input` and `wall_units`
(`m` or `mm`; surrounding spaces and upper-case letters are accepted), and
enable `output_results_netcdf`.

The results then contain `wall_hit_triangle(particle)`: the **zero-based facet
index in STL file order**, with `-1` for markers without an STL hit. Its
`start_index=0`, `stl_file`, and `stl_facet_count` attributes identify the
convention, source mesh, and original number of facets (including zero-area
facets). The global `wall_units` attribute gives the STL length unit. Use the
same STL with its facet order unchanged; no mesh welding, repair, or reordering is
performed. Facets are no longer reoriented: `wall_hit_normal_cart` follows
the STL vertex winding, which need not point outward. The loader checks the
recorded filename (basename, allowing files to be relocated) and facet count.
Older results without `stl_facet_count` retain the basename-only check. No
content hash is recorded; replacing or reordering a mesh with the same name
and count cannot be detected. Always supply the original STL.

Only `wall_hit == 1` contributes to the mesh load. Other losses in `class_lost`
are not assigned to triangles. `zend` supplies the final normalized speed
`p_i = v_i/v0`; thus collisional slowing down contributes its final energy
fraction `e_i = p_i**2` without requiring a different postprocessor.

For nonnegative marker weights `w_i`, total source alpha power `P_alpha`, and
triangle area `A_j` in square meters, the deposited power and flux are

$$
P_j = P_\alpha\frac{\sum_{i:\,\mathrm{hit}(i)=j} w_i p_i^2}
                            {\sum_{i=1}^{N}w_i},\qquad
q_j = P_j/A_j.
$$

The denominator includes **all source markers**, including confined markers
and losses that did not hit the STL. Equal weights give the usual `1/N`
normalization. Supply weights through `marker_weights` (one per source marker)
or a `marker_weight(particle)` NetCDF variable; otherwise all weights are one.
The explicit array takes precedence. These are relative source-power weights,
not absolute watts or weights already multiplied by the final energy. SIMPLE's
current equal-marker output does not write a `marker_weight` variable.

Weight one is power-correct only when births are sampled in proportion to the
physical alpha source, as in the source-weighted births of PR #544
(`startmode = 7`). Other startmodes do not automatically represent that source
with equal weights; supply appropriate source-power weights. Merely using
`1/N` for those modes does not establish a physical wall heat load.

This normalization assumes a common reference birth energy `E0` and
`v0 = sqrt(2 E0/m)`, so `p_i**2 = E_final/E0`. Mixed birth-energy sources need
a separate source-energy convention and are not handled automatically. The
loader reads momentum by the named `phase` dimension, including the native
NetCDF layout `(particle, phase)` seen by Python.

The sum of triangle powers is the deposited fraction of source alpha power;
it equals the full input power only when all source power reaches the wall
without slowing down. No incidence cosine is applied to `q_j`: marker power is
already deposited on the actual triangle area, rather than a projected area.

```python
from pysimple.wall_heat_load import TriangleHeatLoad

load = TriangleHeatLoad.from_netcdf(
    "results.nc", "wall.stl", total_alpha_power_MW=600.0,
)
print(load.summary(threshold_W_m2=1e5))
load.to_vtk("wall_heat.vtk")
load.to_netcdf("wall_heat.nc")
load.to_stl("wall_heat.stl")  # also writes wall_heat.stl.nc
```

`power_W`, `power_density_W_m2`, and `hit_count` have one entry per STL facet.
`power_density_per_MW` is the flux in W/m² for **one MW of source alpha power**:
`1e6 * sum(w_i * p_i**2) / (sum(w_i) * A_j)`. Divide by `1e6` for
MW/m² per MW. All exported geometry is in meters.

Zero-area facets retain their original indices and have zero power and flux
when unhit. They add no area to the metrics. A recorded hit on a zero-area
facet raises an error, even when its marker weight is zero. The CGAL query
tree excludes zero-area facets while preserving the remaining file-order IDs;
the STL export writes zero normals for them.

The metrics use these definitions:

- `peak_W_m2`: maximum triangle flux.
- `p95_W_m2`: area-weighted 95th percentile over the entire mesh, including
  unhit facets. Sort fluxes and select the first whose cumulative area reaches
  95% of total area. This avoids weighting tiny and large facets equally.
- `wetted_area_m2(threshold_W_m2)`: summed area strictly above the threshold;
  the default zero threshold counts positive deposition only.
- `incidence_distribution(bins=18)`: hit counts and deposited watts in angle
  bins spanning `[0, pi/2]` radians. Zero is normal incidence. These are the
  recorded guiding-center chord angles relative to the facet normal, not
  gyrophase-resolved particle velocities.

VTK exports the original facets as cells with `triangle_index`, `area_m2`,
`hit_count`, `power_W`, `power_density_W_m2`, and
`power_density_W_m2_per_MW` arrays for ParaView. NetCDF exports those arrays and
`vertices(triangle, vertex, xyz)` for indexed engineering-tool import. STL has
no standard scalar cell-data representation, so `to_stl` writes geometry plus
a `.stl.nc` sidecar carrying the same facet ordering and fields. A GECKO workflow
can use the STL geometry and join the sidecar fields by `triangle_index`;
automatic import into a particular GECKO version is not provided.

No smoothing is applied. Sparse hits on small triangles can produce large
statistical peaks; inspect `hit_count` and convergence with marker count.
Material response and finite Larmor radius corrections are outside this API.

For a small reproducible example, the CGAL `test_stl_wall_intersection` test
writes `tetra_m.stl` and `tetra_wall_results.nc` in the test build directory.
It records one hit on each of four tetrahedron facets and one no-hit marker.
The Python `test_triangle_heat_load` integration test consumes that fixture
and checks the hit count, triangle indices, hit planes, and deposited power.
