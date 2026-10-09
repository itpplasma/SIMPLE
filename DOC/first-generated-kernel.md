# First Generated Kernel: Splined Field Evaluation

This document records the decision that closes the gap between the
lazy-fortran generation programme and a production SIMPLE physics path
(GitHub issue #515). The lazy-fortran stack generates numerical kernels from
symbolic definitions and emits them to CPU Fortran, OpenMP offload, OpenACC and
native CUDA (tracker: lazy-fortran/fortgen#1). Every issue in that stack ends
inside the toolchain itself; none of them touches code that publishes a physics
number. This page names the first SIMPLE kernel that a generated implementation
must replace, explains why, and fixes the measurement instrument that decides
whether the generated kernel is worth shipping.

The decisive sentence this document sets out to make:

> For this kernel, the arithmetic the mathematics requires is `N_sym`, the
> generated code emits `N_emit`, the compiler produces `N_machine`, the
> hardware delivers `T`, and the accuracy is `X` ulp against a high-precision
> reference.

---

## 1. The choice: splined field evaluation

**The first generated consumer is `splined_evaluate` / `splined_evaluate_with_der`
in `src/field/field_splined.f90`.**

These two subroutines are the leaf of the per-particle orbit-tracing path. Every
guiding-centre step that uses a splined field calls them, they are straight-line
and arithmetic-heavy (no control flow beyond the batch-spline basis sum), and the
derivative variant `splined_evaluate_with_der` has exactly the shape symbolic
automatic differentiation produces well. This makes it the natural first target
for the generation programme.

### Why not the alternatives

1. **Canonical-coordinate transformation kernels** (`field_can_*` family):
   larger, with substantially more control flow and coordinate-system dispatch.
   A poorer first target: the symbolic source is bigger, the emitted code is
   harder to verify, and the win per unit of toolchain effort is smaller.
2. **RK stage evaluation**: already generated inside the toolchain (fortnum),
   but SIMPLE does not consume the generated version, and the RK stages call the
   field anyway, so the field kernel is upstream of them.

### Evidence that this is where the time goes

- The spline path is the hottest evaluation point in the code: it replaces the
  Biot-Savart / VMEC field evaluation on every orbit step.
- libneo#408's 2.16x spline evaluation improvement is direct evidence that this
  path dominates wall time; halving its data (fortnum#84, mixed precision) is
  flagged as the largest available win, consistent with a memory-bound
  straight-line kernel.
- The batch-spline optimization already merged into SIMPLE (v1.4.0, 2.86x
  speedup) measured the same routine.

The choice is a decision recorded here, not an inference from file size: the
baseline instrument below measures `splined_evaluate` on a real NCSX workload,
so the claim "this is the hot path" is testable in-repo.

---

## 2. What "first generated consumer" means

A generated implementation of `splined_evaluate_with_der` must be drop-in
compatible: same argument list, same `BatchSplineData3D` coefficients, same
outputs, evaluated at the same points, and interchangeable with the hand-written
subroutine behind the `splined_field_t` type. The generated kernel competes with
the hand-written one on the four emitted targets (CPU Fortran, OpenMP offload,
OpenACC, native CUDA), and the acceptance gate is:

- **Golden records unchanged.** The numerical gate is the existing one; a
  generated kernel that changes any golden-record output is rejected regardless
  of speed.
- The kernel is verified by readback against its symbolic source
  (lazy-fortran/fortnum#70), checked against the hand-written implementation
  numerically, and measured with the instrument in this document.

---

## 3. The measurement instrument

`test/tests/test_splined_kernel_bench.f90` measures the hand-written kernel on a
realistic NCSX workload (coils field + VMEC reference coordinates). It reports
the two numbers SIMPLE can measure itself:

- **`T`**: wall-clock time per evaluation of `splined_evaluate` (value only) and
  `splined_evaluate_with_der` (value + derivative), best-of-N over repeated
  loops of `n_timing` evaluations.
- **`X` ulp**: maximum ulp deviation of the splined-field `Bmod` from the direct
  Biot-Savart coils field at the same physical point, across a sweep of points.
  The direct coils field is the high-precision reference for this kernel. The
  reported figure is dominated by spline interpolation error (the spline is an
  approximation of the reference field), which is exactly the accuracy the
  kernel delivers as a physics producer.

The other two counts come from the toolchain, not from SIMPLE:

- **`N_sym`**: the arithmetic the mathematics requires, counted in the symbolic
  definition in fortsym.
- **`N_emit`**: the arithmetic the generated code emits, counted by fortsym.
- **`N_machine`**: the machine instructions the compiler produces, counted from
  the generated code (e.g. assembly instruction count).

When a generated kernel lands, the instrument is rerun unchanged and the
generated numbers are inserted into the statement. The statement is the
deliverable: it is the only way to tell whether the generated kernel is faster,
slower, or as accurate as the hand-written one, and that finding belongs in
SIMPLE where the consequence lands.

### Baseline (this repository, gfortran 14, one core)

Measured by `test_splined_kernel_bench.x` on the NCSX coils workload:

```
T (value only)          = ~1.1e-6 s/eval
T (value + derivative)  = ~3.2e-6 s/eval
X ulp (max, Bmod)       = ~2.7e9 ulp vs direct field  (~6e-7 relative)
```

The ulp figure is spline-interpolation accuracy at the default 62x63x64 grid,
not a rounding-error bound. It is the accuracy the kernel delivers today; a
generated kernel must match it (or the golden-record gate rejects it).

---

## 4. Related

- lazy-fortran/fortgen#1 -- tracker for the generation and measurement programme
- lazy-fortran/fortnum#78 -- end-to-end objective `benchmark-orbit-proxima`
- lazy-fortran/fortnum#70 -- readback verification against the symbolic source
- lazy-fortran/fortnum#84 -- mixed precision; the spline path is memory-bound
- lazy-fortran/fortnum#80 -- attribution of the 1.70x CUDA-over-OpenACC gap
- `src/field/field_splined.f90` -- the kernel this document names
- `test/tests/test_splined_kernel_bench.f90` -- the measurement instrument
