# SPECTRE orbits: remaining work

Updated 2026-10-18 against `origin/main` `a041da6`. Delete this file when every
item below is done and verified. History, recovery provenance and the proof
obligations referenced from `concept/*/README.md` are in the
[previous plan](https://github.com/itpplasma/SIMPLE/blob/a041da6/PLAN.md).

## Current state (verified)

- Spline-based, event-resolved full orbit (FO) is the SPECTRE reference; the
  GC/FO hybrid (FO entry at 1.5 radial gyroradii, exit check once per radian)
  is the production path. Cost and loss evidence:
  `concept/spectre_loss_study/results/{boozer-cost,fo-residence,profile}/`.
- libneo #427 (quintic spline kernels) is merged and fetched by default.
- All 24 registered SPECTRE tests pass natively (`make test CONFIG=Fast
  TEST=spectre`, Python tests need a `.venv` with numpy and h5py).
- Issue #528 (no gyrophase root at a current sheet) no longer occurs on the
  production path: entry happens 1.5 gyroradii off the sheet, and
  `test_spectre_fo_interfaces_rk45` on the current-carrying `tok2vol.h5`
  asserts zero FO failures.
- Issue #440 (Level-1 GC interface map) is obsolete: the map was removed in
  #531/#533 (ledger F1, `DOC/spectre-rejected-interface-models.md`) and
  `crossing_level` accepts only 3.

## Remaining items, in order

1. **Midpoint axis-chart momentum (this branch).** After a midpoint endpoint
   flips `(r, theta) -> (|r|, theta + pi)`, `p_theta` was extrapolated from
   old-chart jets (`src/orbit_symplectic.f90`, `attempt_midpoint_step` and the
   multi-volume loop). Accept: `test_axis_crossing_midpoint` fails before the
   fix (relative defect 1.6e-3) and passes after (0), with the full
   non-regression suite and golden records unchanged.
2. **Midpoint Newton failure next to the axis.** The same QA Boozer start
   (`r = 1e-5`, `dtau = 2 pi R/64`) returns `SYMPLECTIC_STEP_MAXITER` one step
   after the crossing at `r = 4.6e-7`; two other starts with `r <= 5e-5` fail
   the same way without a crossing. Accept: trace the non-convergence to its
   source line; the near-axis ensemble of `test_axis_crossing` runs with
   `MIDPOINT` at `dtau = 2 pi R/64` without loss and with bounded energy.
3. **GC energy-shell guard.** Chatter at `npoiper2=256` comes from the absolute
   shell band (5e-4); `npoiper2=512` removes it. Decide between keeping the
   absolute band (then set 512 in SPECTRE examples and tests that trace) and an
   increment test. Accept: 64-marker NCSX hybrid at `npoiper2=256` shows no
   bulk recoveries and unchanged loss within the 95% CI; tok2vol p_zeta
   invariant spread stays below 1e-5.
4. **SPECTRE golden record (#442).** No `test/golden_record` case uses a
   SPECTRE field, and `run_golden_tests.sh` stages only VMEC files. Accept: a
   `tok2vol.h5` hybrid case runs under `make test-golden`, the existing 13
   cases stage unchanged, and the PR states that the case is vacuous until it
   is on `main`.
5. **Poincare validation tolerances (#437).** Accept: `test_spectre_poincare`
   checks against the reference sections at the tolerances in the issue, or
   the issue records measured reasons for looser ones.
6. **FO cost after #427.** Re-profile the hybrid and FO (RK45 stages and GC
   Meiss evaluations), then test larger accepted RK steps at the same
   tolerance. Accept: core-s per particle-s against the `profile` records, with
   matched loss statistics and the tok2vol p_zeta invariant below 1e-5.
7. **Longer and larger validation** at the production default: both fields,
   at least 1,024 markers and longer than 2 ms, with the p_zeta invariant as
   the switch gate. Accept: FO-vs-hybrid final loss fraction within the 95% CI
   and recorded under `concept/spectre_loss_study/results/`.
8. **Close stale issues** #528 and #440 with the evidence above, after
   Chris confirms.

## Deferred, not rejected

Phase-carrying GC switch, second-order GC, encounter maps, resolved GL2,
variational symplectic FO and NCSX 2/4/8-volume equilibria stay in their
`concept/` directories; reopen one only on a measured cost or accuracy gap.
Rejected models F1-F8 stay retired
(`DOC/spectre-rejected-interface-models.md`). Marker 129 and paired NCSX status
differences are diagnostics, not tuning targets.

## Rules

Iterate with runs under a minute (64 markers, 200 us, matched starts); announce
256+ marker 2 ms runs and run them in the background. Use native `make`/CTest.
Do not relax thresholds or regenerate golden data. No protected-host work.
