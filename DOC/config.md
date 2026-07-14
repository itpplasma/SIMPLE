* For symplectic Euler, `npoiper2=256` should be set in particular if all orbits
  (including strongly passing) are traced. This is especially the case when they can collide into this region if `sw_coll=.True.`. If one sets a limit via
  e.g. `contr_pp=-1d0` to trace only mildly passing orbits, or even switch
  them off via `notrace_passing=.True.`, one may resort to 128, or even 64
  in extreme cases with short tracing time. However, also for these cases,
  `npoiper2=256` appears to be a safer choice.

* `boundary_event_fraction_tolerance` sets the dimensionless fractional-step
  bracket tolerance for an `s=1` event. `boundary_event_radial_tolerance` sets
  its dimensionless radial residual tolerance in the integrator coordinate.
  Both default to `-1`, which derives the tolerance from `relerr`. Set positive
  values to refine event location independently of the nonlinear solve.

* `symplectic_newton_warning_mode` defaults to `.true.`. When an implicit solve
  reaches its iteration limit, SIMPLE continues only if the final Newton
  correction is finite and no more than ten times the requested relative
  tolerance. Each occurrence is reported by the corresponding `*_maxit`
  diagnostic. Set it to `.false.` to stop the affected orbit at its last
  accepted position. Larger corrections, exterior-domain states, singular
  linear systems, non-finite values, and unresolved boundary events always
  stop the orbit. Numerical stops are recorded in `orbit_exit_code` and as
  `NaN` in `times_lost`; they are unresolved markers, not physical losses.

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
