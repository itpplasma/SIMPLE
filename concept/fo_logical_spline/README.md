# Superseded W4 spline prototype

The stopped W4 draft was recovered from `SIMPLE-w4` (`45c0598`). Its provider,
generator, benchmark and integration patch were never registered with CMake.
The built and measured replacement is [PR #534](https://github.com/itpplasma/SIMPLE/pull/534),
implemented in `src/spectre_fo_field.f90` and `src/spectre_fo_hybrid.f90`.

The unused draft code and patch were removed during the 2026-10-15 integration
cleanup. They remain in the
[pre-cleanup Git tree](https://github.com/itpplasma/SIMPLE/tree/38ea153/concept/fo_logical_spline),
including the original worker metadata and reported 11-identity generator.
Do not apply that patch over the tested provider. This retires a superseded
implementation, not logical-coordinate splines or the active GC1+/GC2 work.

Current controls, axis-limit evidence, provenance and reproduction commands
are in `concept/spectre_loss_study/` and the root `PLAN.md`.
