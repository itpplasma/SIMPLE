This file contains a (possibly non-exhaustive) list of important changes.

Changes marked with "Y" in the Differing Output column may yield different output when run on the same input
compared to earlier versions.

| Differing behavior | commit hash / tag and timestamp | Summary |
| ---------------- | ------------------------------- | ------- |
| Y | 2026-10-24 | Pin libneo defaults to checked hybrid VMEC Boozer tables and checked EFIT contour conversion for the optional IMAS interface. Fields and orbit results may change, and equilibria that fail quality checks are rejected before tracing. |
| Y | 2026-10-20 | Accelerate the standalone full-orbit Gauss-2 energy/oracle evaluator with exact coefficient sparsity, shared recurrent Fourier phases and FortSym-generated field/Hessian contractions. The algebraically equivalent factored Hessian can change roundoff; orbit length, seeds, step policy, tolerances and energy checks are unchanged. |
| Y | 2026-10-20 | Keep classifier surface Bmin/Bmax per particle. Multi-surface OpenMP classification plots now retain the serial trapped-passing boundary in `fort.20000`; particle classification and loss accounting are unchanged. |
| N | 2026-10-16 | SIMPLE no longer imports the IEEE intrinsic modules; `finite_check` provides `is_finite` and `quiet_nan`. gfortran saved and restored the floating-point environment in every procedure with a local `USE` of a module whose dependencies imported IEEE, including `velo_can` and `integ_to_ref_albert` in the tracing loop. Albert-coordinate tracing is about 1.2x faster; outputs are byte-identical. |
| Y | 2026-10-16 | SPECTRE hybrid tests the full-orbit exit distance once per radian of gyration and at the end of each microstep instead of after every RK45 substep (about 40% of full-orbit time). Exits occur up to 1 rad later; FO-everywhere orbits are unchanged. See `concept/spectre_loss_study/results/profile/`. |
| Y | 2026-10-16 | Default `spectre_fo_entry_larmor` is 1.5 (was 3): SPECTRE multi-volume tracing enters full orbit 1.5 radial gyroradii from an interface or the outer edge and leaves at 2.5. This halves FO time at unchanged loss statistics; the tok2vol p_zeta invariant stays at 8e-8 (entry 1 breaks it). See `concept/spectre_loss_study/results/fo-residence/`. |
| Y | perf/spectre-fo-logical-spline | SPECTRE full orbit uses per-volume geometry/vector-potential splines and logical-position RK45 away from the axis, removing raw spectral evaluation and regular-path stage inversion. `spectre_fo_grid` controls spatial interpolation resolution. Sheet events, Cartesian velocity, GC reconstruction and exit criteria retain their existing meaning; field-grid convergence is required. |
| Y | perf/spectre-fo-logical-spline | Use a regular first-volume axis table for geometry and vector potential. Polar spline derivative error previously produced a wrong axis field and a failed full-orbit recovery with warning-held intervals. Retain cheap spline evaluation and join the potential before taking its curl. |
| Y | 2026-10-14 | Remove the rejected SPECTRE guiding-centre interface maps and their configuration/logging/recovery pipeline. Multi-volume (`Mvol > 1`) orbit configuration now stops with a link to `DOC/spectre-rejected-interface-models.md`, pending the full-orbit interface treatment of PR #530. Single-volume tracing, field/geometry evaluation and start-point generation remain available. Removed code remains at main revision `b113786` (PRs #440, #443, #441). |
| Y | feat/spectre-fo-interfaces | SPECTRE multi-volume tracing uses only `crossing_level = 3` (default) with `integmode > 0`: guiding centres switch to sharp-field RK45 full orbit near interior interfaces. Keep the rejection of fixed-mu levels 0/1 without aliases. Remove the smoothed-layer Boris interface integrator and its integrator/layer/rotation keys; control entry with `spectre_fo_entry_larmor` and RK45 error with `spectre_fo_tolerance`. |
| Y | PR #526 2026-10-13 | Fortran library interface `simple_api`; `simple.x` runs the same phases. `pysimple.init(config="simple.in")` reads simple.in and reproduces `simple.x`; `init` now starts every run from the compiled defaults and rejects unknown keys. Python calls without `integrator` use the configured `integmode` (previously always midpoint). Removed the unused Python batch layer (`ParticleBatch`, `BatchResults`, samplers). Old input files with the removed keys `local`, `loopskip`, `bmod_ref` or `macrostep_log_ratio` are read again; the keys are ignored with a warning. `generate_start_only` now also stops TEST-field runs. `examples/simple_full.in` lists every key at its default. |
| Y | PR #527 2026-10-13 | Fix the sign of the curvature drift (`hcurl`) for SPECTRE fields and for unsplined reference fields on the RK45 path: both passed component-first `d h_i/dx_j` to `compute_hcurl`, which takes the derivative index first. A forbidden Level-1 crossing that relocates to a transmission now enters the target volume instead of staying owned by the home volume. The tok2vol test fixture is replaced by a valid current-carrying equilibrium; the current-free sheet fixture tok2vol_sheet serves the pure-magnitude interface tests. |
| Y | PR #503 2026-07-22 | In warning mode, recover a failed conventional guiding-centre interval with toroidal regularization and, only at an unresolved polar axis, a Cartesian full-orbit bridge. Preserve the last accepted momentum shell, return to the requested symplectic method after a safe streak, and retain the bounded isolated-warning hold instead of freezing a marker indefinitely. Successful legacy steps are unchanged. |
| N | v1.5.1 (release/1.5) 2025-12-23 | Use atomic updates when incrementing the confined pass/trap counters so each branch stays thread-safe without enclosing the entire branch in an OpenMP critical region. |
| N | v1.4.2 (release/1.4) 2025-12-23 | Backport the atomic-update behavior for confined counters so release/1.4 matches v1.5.1’s threading sync. |
| Y | 6ae316f 2025-12-22 | Deterministic mode: change RNG seed from 0 to 42 |
| Y | 7f24b57 2024-12-06 | Change hcov computation order from matmul(B,J)/|B| to matmul(B/|B|,J) - ~1e-12 rounding differences in coils field |
| N | 1c38cfa 2024-10-16 | Buffered NetCDF orbit output system - replaces fort.9000X files with orbits.nc |
| N | v1.4.0 (493afe3) 2024-09-05 -> `release/1.4` | Batch spline optimization (2.86x speedup), Meiss-Hazeltine and Albert coordinates |
| Y | 5de27cd 2024-08-04 | Change sign of VMEC toroidal flux (v1.3.1 on release/1.3: ac3ba6b)
| Y | e25ca7e 2024-08-24 | Remove zzg and rng.f - starting pitch gets double precision random number xi now |
| Y | c2bb751 2024-08-19 | Restore behavior of `release/1.3` branch |
| Y | 6c5eb2c 2024-03-11 | Different `bmin` and `bmax` from local init lead to 1e-4 difference in trap_par |
| Y | 081a859 2023-12-15 | Set default values to -1d10 for end angles to mark it's not VMEC coordinates |
| Y | 8e19626 2022-06-29 -> `release/1.3` | Removed bmod_ref that rescaled to match B_00 on flux surface |
| N | b5fc202 2022-06-29 | MPI working again, misses later columns in times_lost.dat |
| N | c796d90 2022-06-28 | Merge collisions branch back to master, cross-checked results |
| Y | 7ca8fad 2022-06-23 | Default to old VMEC near-axis healing again |
| N | 86c3a5d 2022-06-23 | simple.in is now in namelist format |
| Y | e264ef7 2022-05-31 | New healing procedure for VMEC near axis |
| Y | 468e317 2022-02-14 | Compute actual momentum in z(4) every timestep |
| Y | a247433 2022-02-13 | Disable extrapolation entirely in init_sympl |
| N | 6981791 2021-12-20 | start.dat is now always stored in VMEC coordinates. |
| Y | v1.1 2020-02-13 | Lambda is now evaluated on half mesh, improves B-field interpolation. |
