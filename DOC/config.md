# Configuration notes

Guidance for choosing `&config` values. Every key, its default and its meaning
are listed in `examples/simple_full.in`.

* For symplectic Euler, `npoiper2=256` should be set in particular if all orbits
  (including strongly passing) are traced. This is especially the case when they can collide into this region if `swcoll=.True.`. If one sets a limit via
  e.g. `contr_pp=-1d0` to trace only mildly passing orbits, or even switch
  them off via `notrace_passing=1`, one may resort to 128, or even 64
  in extreme cases with short tracing time. However, also for these cases,
  `npoiper2=256` appears to be a safer choice.

* `boundary_event_fraction_tolerance` sets the dimensionless fractional-step
  bracket tolerance for an `s=1` event. `boundary_event_radial_tolerance` sets
  its dimensionless radial residual tolerance in the integrator coordinate.
  Both default to `-1`, which derives the tolerance from `relerr`. Set positive
  values to refine event location independently of the nonlinear solve.

* `symplectic_newton_warning_mode` defaults to `.true.`. When an implicit solve
  reaches its iteration limit, the generic integrators accept a finite final
  correction through 100 requested-tolerance units; SPECTRE keeps a 10-unit
  bound. Midpoint steps that cross the polar-coordinate axis have one additional
  local rule. If both radial Newton variables remain within `1d-6`, the endpoint
  radius is negative, and every angular and momentum correction satisfies the
  usual warning bound, SIMPLE accepts the step and applies the equivalent
  positive-radius chart switch. This avoids rejecting a physically regular axis
  passage because its signed radial coordinate is singular. The
  `warning_axis_crossing_accept` diagnostic counts these steps.

  Other rejected steps use bounded dyadic symplectic retries and then the
  established adaptive-RK pusher from the last accepted state. The RK fallback
  is capped at 10,000 internal steps. If both methods fail, warning mode rolls
  back and consumes one isolated microstep, records `warning_step_skip`, and
  retries. A successful step resets the allowance; a second consecutive failure
  is a numerical exit. Set the option to `.false.` for strict diagnostic runs
  that end the affected marker at the first exhausted recovery. Numerical exits
  use codes 101--105 and `NaN` in `times_lost`.

  A collisionless guiding-centre marker that exhausts every recovery method
  after the warning hold while its last validated normalized toroidal flux is
  below `s=0.01` is classified as numerically confined. For chartmap inputs,
  SIMPLE converts the stored `rho` to `s=rho^2` before applying the gate. It
  uses exit code 4 and `trace_time` in `times_lost`. This core-only fallback is
  disabled for collisions, strict mode, full-orbit tracing, SPECTRE, and
  failures outside the axis-local region.

* `canonical_grid_nr`, `canonical_grid_ntheta`, and `canonical_grid_nphi`
  control the Meiss or Albert canonical-map grid. Their defaults are 62, 63,
  and 64. Each dimension must be between 6 and 512, and their product must not
  exceed 2,097,152 grid points. This limit guards allocation arithmetic; the
  requested host and accelerator memory must still fit the selected grid and
  spline backend. `canonical_ode_relerr` controls the radial transformation
  solve and defaults to `1d-11`. NVHPC
  builds use at least `1d-8` only on slices where the covariant toroidal field
  component is near zero, preventing step-size underflow at that coordinate
  singularity. Change one setting at a time in resolution studies.

* Deterministic collisions seed each particle by its original one-based source
  row, rather than repeating one kick sequence for every particle. Untagged
  initialization calls retain the legacy `ran_seed` state. `startmode=4` uses
  the selected `start.dat` row numbers persisted as text in `batch.dat`.
  Only mode 4 applies that selection; other modes retain their existing
  `batch_size` truncation behavior. Batch draws use `ran_seed` even when
  deterministic collision replay is disabled.
  For a contiguous external shard, set `particle_index_offset` to the number
  of original source rows preceding that shard (default zero). The same marker
  then replays the same stream across thread counts, native selection and
  contiguous shards. Reordered or noncontiguous external subsets must retain
  their original ordering/identity rather than treating local row numbers as
  global IDs. Offset plus source row must be between 1 and 2147483646.
  Different deterministic streams are not a proof of statistical independence.

* The full Coulomb operator subcycles each collision interval using coefficients
  at the current momentum. `stost`'s optional `max_substep_fraction` (default
  0.1, supported range 0 < fraction <= 0.25) controls numerical resolution:
  `dhh*h <= fraction^2`, `abs(fpeff)*h <= fraction*p`, and, when energy
  diffusion is active, `sqrt(6*dpp*h) <= fraction*p`. These bounds are selected
  before random draws. They introduce no thermalization or energy cutoff.
  Pitch diffusion uses an isotropic random rotation through
  `sqrt(4*dhh*h)` rather than an Euler pitch kick followed by reflection or a
  random reset. Its generator is the same Lorentz pitch operator, to first
  weak order in the substep. Energy diffusion and drift retain their existing
  generator and uniform, unit-variance radial noise. A substep underflow or
  the 100000-substep budget limit returns an explicit error with the original
  state restored. Exact vacuum leaves both the state and RNG untouched.
