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
- Midpoint axis chart switch (this branch, `fix/spectre`): after
  `(r, theta) -> (|r|, theta + pi)` the field is evaluated fresh instead of
  extrapolating `p_theta` from old-chart jets; `test_axis_crossing_midpoint`
  defect 1.6e-3 before, 0 after; 150/150 non-regression tests pass. The
  crossing that this test produces is a warning-mode acceptance, not a
  converged root (item 1).

## Remaining items, in order

1. **Midpoint step across the axis: needs a regular near-axis chart
   (decision for Chris).** Traced on `fix/spectre` 047da79. The probes and
   the near-axis ensemble test are kept outside the repository with the
   handoff.
   - Verified: every failing solve has the same 2-cycle. The failing solves
     are the `test_axis_crossing_midpoint` start, identical to 1%-radius
     `test_axis_crossing` start 3, at step 156, and scaled starts 4 and 6 at
     steps 1255 and 2. The first Newton update lands at `r < 0`, the floor
     `if(x(1) < 0.0) x(1) = 0.01` (`src/orbit_symplectic.f90:827-828`) moves
     the iterate to `r = 0.01`, where the radial residual is 1.8e5, and the
     next update lands at `r < 0` again. Angles and `p_phi` converge for
     that floored radius; `r` never does.
   - Verified: the accepted crossing in `test_axis_crossing_midpoint` is not
     a root of the midpoint equations. `accept_warning_axis_crossing`
     (`:89-110`, called at `:875`) checks only components 2-4 and accepts the
     update from the floored iterate (`xlast(1) = 0.01`, final
     `x(1) = -4.64e-7`). The next step then fails. The #503 unit test
     `test_newton_solver_status` encodes exactly this floored state as an
     accepted crossing.
   - Verified (multi-start numerics, not a proof): the failing steps have no
     root with `r_end, r_mid > 1e-11`. Seeds tried: 3072 for the first start,
     covering all theta_end in pi/32 steps, `r_end` from 1e-8 to 8e-5 and
     `r_mid` from 1e-8 to 4e-5, with positivity-damped Newton; 3072 each for
     the other two failing steps. Euler1 fails at the same three steps.
     Retry subdivision to depth 8 still fails on start 6 and accepts 41
     warning-mode crossings, with `max|dH|/H` up to 6e-4.
   - Hypothesis (mechanism): in the `(s, theta)` chart,
     `dB/ds = (dB/drho)/(2 rho)` (libneo `boozer_converter.F90:337`; it
     evaluates at `abs(r)`, `:304`). The midpoint theta residual therefore becomes singular as `r_mid -> 0`, and
     a step whose chord passes the axis closer than about one step length
     has no discrete solution. This is a property of the chart, not of
     Newton conditioning or tolerances.
   - Rejected: damping the step so that `r` stays positive. It drives `r` to
     about 1e-20, where the absolute radial tolerance (`tolref(1) = 1`)
     accepts collapsed non-roots.
   - Options: (a) take the symplectic step near the axis in a regular
     canonical chart (`X, Y = sqrt(2 p_theta) (cos, sin) theta`), a new
     scheme; or (b) reject the floor-cycle state, so that the step returns
     MAXITER and the production RK axis recovery (`orbit_timestep_axis`,
     regular chart) takes over. Option (b) reverses the #503 acceptance and
     can change near-axis golden markers.
   - Accept: with the chosen option, the 1%-radius ensemble runs with
     `MIDPOINT` at `dtau = 2 pi R/64` with no loss or fabricated crossing and
     with bounded energy, and the golden records stay unchanged.
2. **GC energy-shell guard.** Chatter at `npoiper2=256` comes from the absolute
   shell band (5e-4); `npoiper2=512` removes it. Decide between keeping the
   absolute band (then set 512 in SPECTRE examples and tests that trace) and an
   increment test. Accept: 64-marker NCSX hybrid at `npoiper2=256` shows no
   bulk recoveries and unchanged loss within the 95% CI; tok2vol p_zeta
   invariant spread stays below 1e-5.
3. **SPECTRE golden record (#442).** No `test/golden_record` case uses a
   SPECTRE field, and `run_golden_tests.sh` stages only VMEC files. Accept: a
   `tok2vol.h5` hybrid case runs under `make test-golden`, the existing 13
   cases stage unchanged, and the PR states that the case is vacuous until it
   is on `main`.
4. **Poincare validation tolerances (#437).** Accept: `test_spectre_poincare`
   checks against the reference sections at the tolerances in the issue, or
   the issue records measured reasons for looser ones.
5. **FO cost after #427.** Re-profile the hybrid and FO (RK45 stages and GC
   Meiss evaluations), then test larger accepted RK steps at the same
   tolerance. Accept: core-s per particle-s against the `profile` records, with
   matched loss statistics and the tok2vol p_zeta invariant below 1e-5.
6. **Longer and larger validation** at the production default: both fields,
   at least 1,024 markers and longer than 2 ms, with the p_zeta invariant as
   the switch gate. Accept: FO-vs-hybrid final loss fraction within the 95% CI
   and recorded under `concept/spectre_loss_study/results/`.
7. **Close stale issues** #528 and #440 with the evidence above, after
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
