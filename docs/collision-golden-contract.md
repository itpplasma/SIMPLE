# Collision regression after independent marker streams

The `collision` golden fixture retains its original 32 markers, Boozer field,
nonzero flat collision coefficients, deterministic starts and 0.1 ms horizon.
Its run and the ordinary historical runtime comparison remain required.
Only this case explicitly selects `compare_collision_results.py` through
`GOLDEN_RECORD_STOCHASTIC_COLLISION=1` in CTest. All other cases and direct
comparisons retain their existing strict defaults and tolerances.

Previously, resetting the same seed at each marker correlated the collision
kicks. Correcting marker stream identities intentionally changes individual
stochastic histories. The corrected bounded angular diffusion and low-momentum
substeps can also change outcomes. Particlewise historical equality is therefore
not the acceptance criterion for this one fixture. The old 32 correlated
histories are not an independent physical distribution oracle.

The replacement checks all 32 marker IDs, unchanged starts and initial
invariants, finite current outcomes, physically admissible endpoints, positive
momenta, bounded pitch, matching time/status files and correct trace-horizon
accounting. Only completed and LCFS exit codes (0/1) are allowed. WALL code2
is impossible for this no-wall fixture and fails, as does a nonblank
`wall_input` in either reference or current input. The compiled `wall_input`
default is empty; `wall_units='m'` does not enable a wall. Numerical failures,
skipped markers and core classifications also fail.
With at least two survivors, roundoff-identical momenta fail the repeated-kick
regression check for this nonzero-coefficient, equal-start-momentum fixture.
The comparator prints changes in historical loss count and mean resolved p²
as descriptive evidence; no acceptance bound is calibrated against the old
correlated sample. This limited check is not independent distribution validation.

Scientific acceptance requires the independent checks already registered in
the full native test/CI authorities:

- `test_collision_substeps`: pitch P1/P2 moments, heat-kernel convergence,
  radial Ornstein–Uhlenbeck moments, Maxwellian equilibrium and exact vacuum;
- `test_particle_rng`: distinct marker streams, seeded replay and unchanged
  streams across OpenMP scheduling, native selection and external shards;
- `test_collisional_stl_wall`: the open CGAL pipeline fixture, including seeded
  repeatability and recorded slowed-down wall-hit energies.

Passing the limited golden fixture alone is insufficient. The complete PR
must retain these physics checks, the other golden cases and the performance
gate. It must record the measured collision loss/energy movement and obtain
reviewer acceptance under `docs/review-policy.md`. There is no predetermined
sign for movement caused by stream decorrelation. The low-momentum correction
must be explained separately using its analytic weak-moment/convergence evidence.
No reference patch, global tolerance change or skipped-case exception is used.
