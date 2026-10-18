# Numerical settings

These `&config` keys set step size, event location, failure handling and
random streams. Change them only for resolution studies or diagnostics; the
defaults in `examples/simple_full.in` suit loss-fraction runs. Change one
setting at a time and compare loss statistics, not single trajectories.

## Step size

| Key | Default | Meaning |
|---|---|---|
| `npoiper2` | 256 | microsteps per path length $2\pi R_0$: step $\Delta\tau \le 2\pi R_0/N$ (cm) with $N$ = `npoiper2`, reduced so that each macrostep holds a whole number of steps |
| `relerr` | `1d-13` | relative tolerance of the Newton solves and of RK45 |

$R_0$ is the major radius of the equilibrium (`src/params.f90`, `params_init`).
For symplectic Euler, keep `npoiper2 = 256` when strongly passing orbits are
traced, and always with collisions (`swcoll = .true.`), which scatter markers
into the passing region. Without collisions and with strongly passing orbits
excluded (`contr_pp`, `notrace_passing`), 128 is usually enough, 64 only for
short `trace_time`.

To check a coarser step, run the same `start.dat` at two values of `npoiper2`
and compare the loss statistics at the full trace time:

```bash
python tools/step_convergence.py coarse/times_lost.dat fine/times_lost.dat
```

The coarse step is adequate when the confined-fraction bias and the loss-time
distance are below the Monte Carlo error, about $1/\sqrt{\text{markers}}$.

## LCFS event location

A symplectic step that leaves $s = 1$ is shortened to the crossing by a
bracketed root search.

| Key | Default | Meaning |
|---|---|---|
| `boundary_event_fraction_tolerance` | `-1d0` | width of the bracket in units of one step; `-1` uses `max(1d-12, 10*relerr)` |
| `boundary_event_radial_tolerance` | `-1d0` | radial residual in integration coordinates; `-1` uses `max(1d-10, 10*relerr)` |

Positive values refine the event independently of `relerr`
(`src/orbit_symplectic.f90`).

## Failed implicit steps

`symplectic_newton_warning_mode` (default `.true.`) decides what happens when
a Newton solve reaches its iteration limit.

- `.true.`: SIMPLE accepts a finite final correction within 100 times the
  requested tolerance (10 times for SPECTRE). Near the polar axis it also
  accepts a step whose endpoint radius changes sign when both radial Newton
  variables stay below $10^{-6}$ and all angular and momentum corrections meet
  the same bound; the step continues in the equivalent positive-radius chart
  (counter `warning_axis_crossing_accept`). A step that still fails goes
  through the recovery sequence in {doc}`/model/integrators`. A marker whose
  recovery is exhausted ends with a numerical exit code (101–105, see
  {doc}`/user-guide/outputs`).
- `.false.`: the marker ends at the first exhausted recovery. Use this for
  strict diagnostic runs.

In warning mode, a collisionless guiding-centre marker of a symplectic run
whose recovery is exhausted while its last valid $s < 0.01$ counts as
confined, with exit code 4 and `trace_time` as loss time
(`src/simple_main.f90`, `should_classify_numerically_confined`). Chartmap runs
convert $\rho$ to $s = \rho^2$ first. This rule is off with collisions, in
strict mode, for full orbits and for SPECTRE.

## Canonical-map grid (Meiss, Albert)

| Key | Default | Meaning |
|---|---|---|
| `canonical_grid_nr`, `canonical_grid_ntheta`, `canonical_grid_nphi` | 62, 63, 64 | grid of the canonical transformation and field splines; each 6–512, product $\le 2\,097\,152$ |
| `canonical_ode_relerr` | `1d-11` | tolerance of the radial transformation ODE |

The product limit guards allocation size; host and device memory must still
hold the quintic splines. NVHPC builds raise the ODE tolerance to at least
$10^{-8}$ on slices where the covariant toroidal field component is near zero,
to avoid step-size underflow there (`src/field/field_can_meiss.f90`).

## Random streams

With `deterministic = .true.`, each marker draws from its own stream, seeded
from `ran_seed` and the marker identity $n = n_0 + i$: $n_0$ is
`particle_index_offset` and $i$ the one-based row in `start.dat` (for
`startmode = 4`, the original row stored in `batch.dat`). The same marker
receives the same collision kicks for any thread
count and for native or external sharding. For a contiguous external shard,
set `particle_index_offset` to the number of rows before the shard (default
0); $n$ must lie in $[1, 2\,147\,483\,646]$. A reordered or non-contiguous
subset must keep the original row numbers. Distinct streams are not a proof
of statistical independence.

Batch draws (`batch_size`) use `ran_seed` also without `deterministic`.
Collision substeps are chosen automatically; their bounds are in
{doc}`/model/collisions`.

## Particle energy, Lorentz collisions and radial potential

- `particle_energy_eV` (default `-1`) sets the energy directly. With `-1` the
  energy is `3.5d6/facE_al` eV. A nondefault `facE_al` must agree with an
  explicit energy, otherwise input validation fails.
- `collision_model` (default `'full'`): `'lorentz_nustar'` selects
  monoenergetic pitch-angle scattering for `swcoll = .true.`; requires
  `nu_star_standard >= 0`, `lorentz_major_radius_cm > 0` and nonzero
  `lorentz_iota`, with `nu_star_standard = R0*nu/(abs(iota)*v0)`. `p_abs`
  stays unchanged. See {doc}`/model/collisions`.
- `radial_electric_potential_slope` (default `0`) adds `Phi(s) = slope*s` to
  the Boozer Hamiltonian.
- With `output_results_netcdf` and `output_orbits_macrostep`, `results.nc` and
  `orbits.nc` carry the effective energy, collision model and these values as
  global attributes.
