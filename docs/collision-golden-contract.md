# Stochastic collision golden case

The golden case `test/golden_record/collision` traces 32 markers with
collisions in Boozer coordinates for 0.1 ms, with the flat default `&profiles`
and deterministic starts. Each marker has its own random stream, so individual
histories are not reproducible across changes of the collision operator or of
the stream assignment, and the case cannot be compared particle by particle.
CTest therefore sets `GOLDEN_RECORD_STOCHASTIC_COLLISION=1` for this case
only, which selects `test/golden_record/compare_collision_results.py`. All
other cases keep the strict comparison.

The comparator accepts the current run when:

- `simple.in` sets `ntestpart = 32`, `swcoll = .true.`,
  `deterministic = .true.`, no `wall_input`, and the same trace time as the
  reference;
- all 32 marker IDs are present in order, starts and initial invariants equal
  the reference, and `times_lost.dat` agrees with `orbit_exit_code.dat`;
- every exit code is 0 (completed) or 1 (LCFS); wall, skipped, core and
  numerical exits fail;
- all outcomes are finite, momenta positive, pitches in $[-1, 1]$, LCFS
  endpoints near $s = 1$, and loss times inside the trace horizon;
- no two surviving markers have momenta equal to round-off, which would
  indicate a repeated random stream.

It prints the change of loss count and mean $p^2$ as information; there is no
bound on it. This check is not a test of the collision physics. That is
covered by:

- `test_collision_substeps`: pitch moments, heat-kernel convergence, radial
  Ornstein–Uhlenbeck moments, Maxwellian equilibrium, exact vacuum;
- `test_particle_rng`: distinct marker streams, seeded replay, identical
  streams across OpenMP scheduling and external shards;
- `test_collisional_stl_wall`: collisional tracing to a CGAL wall, seeded
  repeatability and slowed-down wall-hit energies.

A pull request that changes collision results must pass these tests and state
the measured change of losses and energy (see the contributor guide in
{doc}`/contributing`).
