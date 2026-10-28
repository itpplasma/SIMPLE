# Heat-load estimators

This page defines how SIMPLE's post-processing turns lost markers into
deposited power and heat flux, on flux-surface bins or on STL wall triangles,
and what limits the accuracy of these estimates.
Usage: {doc}`/user-guide/wall-heat-loads`.

## Marker power and energy

For D-T, $P_\alpha = P_\text{fus}/5$. Marker $i$ represents the share
$w_i/\sum_{k=1}^N w_k$ of the source alpha power, with relative weights
$w_i \ge 0$; equal weights give $P_\alpha/N$. Equal weights represent the
physical source only when births are sampled in proportion to it
(`startmode = 7`, {doc}`/alpha-birth`).

All markers start with the reference energy $E_0 = 3.5\,\text{MeV}/$`facE_al`,
$v_0=\sqrt{2E_0/m}$. At the exit the marker carries

$$
\frac{E_i}{E_0} = p_i^2, \qquad p_i = \frac{v_i}{v_0} = \texttt{zend(4)} .
$$

Without collisions $p_i = 1$. With collisions the energy given to the
background before the exit is not deposited on the wall; $p_i < 1$ on average,
but energy diffusion can raise individual markers above $p_i = 1$.

Mixed birth energies (for example the full, half and third energy components
of a neutral beam) start with $p_{0,i} = $ `zstart(4)` $\neq 1$. The STL facet
estimator below then treats $w_i$ as relative numbers of source particles and
normalizes by the source energy, $\sum_k w_k p_{0,k}^2$, so that a
collisionless run in which every marker reaches the wall deposits exactly the
source power. For monoenergetic sources ($p_{0,k}=1$) this reduces to
$\sum_k w_k$. The flux-surface bin estimator still assumes a monoenergetic
source.

The particle and energy loss fractions and the lost power are

$$
f_\text{particle} = \frac{n_\text{lost}}{N}, \qquad
f_\text{energy} = \frac{\sum_{i\in\text{lost}} w_i p_i^2}{\sum_{k} w_k}, \qquad
P_\text{lost} = f_\text{energy} P_\alpha .
$$

## Flux on flux-surface bins

Lost markers that crossed the LCFS (exit code 1) are binned by their exit
angles $(\vartheta_i, \varphi_i)$ = `zend(2:3)`. With equal weights

$$
q_\text{bin} = \frac{P_\alpha}{N A_\text{bin}} \sum_{i \in \text{bin}} p_i^2 .
$$

The bin area comes from the outermost chartmap surface
$\mathbf r(\vartheta,\zeta)$ (the LCFS for a chartmap written by
`vmec_to_chartmap.x`), stored for one field period and repeated $n_\text{fp}$
times:

$$
A_\text{bin} = \iint_\text{bin}
\left|\frac{\partial \mathbf r}{\partial \vartheta}\times
\frac{\partial \mathbf r}{\partial \zeta}\right| \mathrm d\vartheta\,\mathrm d\zeta ,
\qquad
A = n_\text{fp}\int_0^{2\pi}\!\!\int_0^{2\pi/n_\text{fp}}
\left|\partial_\vartheta \mathbf r\times\partial_\zeta \mathbf r\right|
\mathrm d\zeta\,\mathrm d\vartheta .
$$

The integral is a sum over chartmap grid cells whose nodes lie in the bin, so
a bin needs grid points: use at most `ntheta` poloidal and
`nzeta`$\times n_\text{fp}$ toroidal bins of the chartmap, otherwise bins get
zero area. Bin angles and chartmap angles must be the same coordinates (VMEC
angles for `vmec_to_chartmap.x`).
Implementation: `compute_bin_areas_from_chartmap` in `python/pysimple/engineering.py`.

## Load on STL triangles

A marker with `wall_hit == 1` deposits on its facet $j$ = `wall_hit_triangle`.
For facet area $A_j$ (m²)

$$
P_j = P_\alpha\,\frac{\sum_{i:\,\text{hit}(i)=j} w_i\, p_i^2}{\sum_{k=1}^{N} w_k\, p_{0,k}^2},
\qquad q_j = \frac{P_j}{A_j}.
$$

The denominator runs over **all** source markers, including confined markers
and losses that did not reach the STL; $p_{0,k}$ = `zstart(4)` is one when
`results.nc` has no `zstart`. $\sum_j P_j$ is the deposited part of the source
power $P_\alpha$. No incidence cosine is applied: the power lands on the actual
facet area, not a projected one. Weights come from `marker_weights`, else
from a `marker_weight(particle)` variable in `results.nc` (SIMPLE does not
write one), else all ones. Facets of zero area keep their index, receive no
power and add no area; a recorded hit on one is an error.
Implementation: `TriangleHeatLoad.from_netcdf` in `python/pysimple/wall_heat_load.py`.

Metrics of `TriangleHeatLoad.summary`:

| Metric | Definition |
|---|---|
| `peak_W_m2` | $\max_j q_j$ |
| `p95_W_m2` | area-weighted 95th percentile over all facets, unhit included: sort $q_j$, take the first whose cumulative area reaches 95 % of the total |
| `wetted_area_m2(t)` | $\sum A_j$ over facets with $q_j > t$ (default $t = 0$) |
| `incidence_distribution(bins=18)` | hit count and deposited power in angle bins on $[0, \pi/2]$ |

The incidence angle is $\arccos|\hat{\mathbf s}\cdot\hat{\mathbf n}|$ between
the guiding-centre step chord $\hat{\mathbf s}$ and the facet normal
$\hat{\mathbf n}$ (0 is normal incidence), not a gyrophase-resolved velocity.
The normal follows the STL vertex winding and need not point outward. The
loss time is interpolated linearly along the intersecting chord
(`locate_wall_segment` in `src/simple_main.f90`).

## Accuracy

- **Counting noise.** A bin or facet with $k$ hits has relative error about
  $1/\sqrt{k}$. No smoothing is applied; small facets with few hits produce
  large peaks. Converge peak values in marker number.
- **Guiding centre, not impact point.** Positions are guiding-centre
  positions. The particle hits up to one gyroradius away,
  $\rho_\alpha = m_\alpha c\, v_\perp/(Z_\alpha e B) \le 5.4$ cm for 3.5 MeV
  alphas at 5 T. Fine wall features and peak locations below this scale are
  not resolved.
- **Hit energy.** With collisions, the hit keeps the momentum of the kicks
  applied before the intersecting microstep; the energy has the time
  resolution of the microstep.

## References

- C. G. Albert et al., Symplectic methods for optimization of wall heat loads
  by fusion alpha particles in stellarator reactors, Nucl. Mater. Energy 48,
  102200 (2026), [doi:10.1016/j.nme.2026.102200](https://doi.org/10.1016/j.nme.2026.102200).


### Structured event output

`results.nc` schema 1.1 adds `kinetic_energy_eV` for every final marker state,
using the run's effective reference energy times `zend` momentum squared.
For first wall events select both `wall_hit == 1` and `orbit_exit_code == 2`;
`times_lost` is the event time in seconds. Other final states are not impacts; untraced markers retain the sentinel state
from `zend`, so their stored energy is not a physical impact energy.
The existing position variables carry their length units; Cartesian positions
and facet normals share the right-handed field/chart and STL mesh frame.
Facet indices are zero-based in STL file order, with -1 for no hit. Particle
arrays retain source order, without claiming a globally unique marker ID.

`marker_weight` contains dimensionless relative source weights. By default all
weights are one and its `source` attribute explicitly says that sampling weights
are unknown. This default is appropriate for equal-marker summaries, and does
not establish a physical NBI source normalization. Callers of the Fortran
`write_results_netcdf` routine may supply `marker_weights`; post-processors may
instead supply the original external weights. Normalize over the complete
source population, including confined markers, and supply physical source power
or particle rate separately.

Pitch remains `v_parallel/v`, with speed normalized by the global `v0`.
Incidence is the absolute angle between the guiding-centre integration chord
and the wall facet normal. Neither gyrophase nor a physical impact velocity
vector is recorded, and pitch does not determine that vector. Finite Larmor
radius impact displacement and instrument/thermal response are not represented.

The output stores background ion masses and charges, the complete configured
`profiles` namelist and `profile_label_file`. When a tabulated profile is active,
`collision_profile_table` stores its actual samples, including tables supplied
through the library API. Temperatures are in eV and densities in m^-3. Analytic
profiles use the scales and coefficients in the namelist. External field,
chart, wall and profile-label files still need to be archived with the run;
filenames are not content identities. Numeric arrays use lossless NetCDF4
shuffle/deflate compression, without changing their floating-point precision.
