# PR #530 rebase regression repair

- Starting branch: feat/spectre-fo-interfaces, 9f5f395; working tree clean.
- Read CLAUDE.md and relevant Fortran guidance; using make only, no delegation.
- review_h.md is absent; supplied findings checked against production source.
- Confirmed simple.x and simple_api share init_spectre_field.
- Existing direct volume test overrides sympl_rmax and masks the init regression.
- Added production-init volume-2 GC witness without overriding sympl_rmax.
- Fail before (9f5f395 sources): bound 1; first step enters FO; angular advance 0.
- Pass after: bound 3; 10000 steps; FO entries/stops/warning holds all zero;
  angular advance 47.010 rad. Outward edge marker loses at rho_g=2,
  status 88, terminal fraction 3.1826e-4; one landing, no FO or stops.
- Restored Mvol+1 trial margin; reset default bound on field initialization.
- Updated bound comments and coordinates documentation. Interface doc has no
  bound claim and needs no edit.
- Fix committed as 32e5c1e.
- Removed unread stop-reason counters and unused marker/time argument chain;
  warning-hold reason/interface/direction matching retained unchanged.
- Correctness contract: T3 init fix; rho_g is the dimensionless global volume
  label, angles radians. Physical edge and event localization remain at Mvol.
  Cleanup is T1; counters still consumed by diagnostics retain their behavior.
- Regression uses actual GC motion and absence of recovery, not bound equality.
  Edge witness checks physical loss and fractional terminal time. No tolerances
  or golden records changed.
- make CONFIG=Fast passed before and after the changes.
- Full make test CONFIG=Fast VERBOSE=1: 95/95 passed, no failures or skips,
  765.02 seconds. Includes FO interface/events, invariant, loss, Python API,
  golden-orbit output and convergence tests; no tolerances loosened.
- Full-suite new witness after cleanup matches pre-cleanup values exactly:
  10000 steps, 0 FO entries/stops/holds, 47.010 rad; loss at 2.0000000000,
  status 88, fraction 3.1826e-4. Logs: build/h2_fail_before.log,
  build/h2_pass_after.log, build/h2_full_test.log.
- Inspected diff; only requested fix, cleanup, focused test, docs and this record.
- Remaining compiler real-comparison warnings are on unchanged source lines.
- Publication target: plain push to origin/feat/spectre-fo-interfaces;
  no merge, GitHub post, force push, other agent, or protected-host work.
