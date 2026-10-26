# Coding Standard

This standard records the conventions used in SIMPLE. Rules marked
**(new code)** are not yet met everywhere in the existing sources; they apply
to new code and to code that a pull request rewrites. Do not reformat or
restructure unrelated code in the same pull request: a formatting-only or
modernisation change is its own pull request (tier T1, golden records
unchanged).

## Fortran

### Language and files

- Modern free-form Fortran with lower-case keywords, compiled with GNU
  Fortran. Code on the OpenACC paths must also build with NVIDIA Fortran
  (`make nvfortran`, `make nvfortran-acc`). Native CUDA requires the CUDA
  toolkit and either GNU or NVIDIA Fortran; see
  [the GPU build instructions](../user-guide/gpu.rst).
- All Fortran sources pass through the preprocessor (`-cpp`), whatever the
  file extension.
- The GNU build uses `-Wall -Wextra -Wimplicit-interface` and rejects lines
  longer than 132 characters. Existing code emits warnings; a change must not
  add new ones.
- One module per file, with the file named after the module **(new code)**.

### Modules and declarations

- Every procedure lives in a module; no new external procedures or `include`
  files **(new code)**.
- Each module starts with `use` statements, then `implicit none`, then
  `private` with an explicit `public` list for new modules.
- Use `use module, only: name` for every module dependency **(new code)**.
  Rename on import to avoid clashes, for example
  `use field_can_mod, only: eval_field => evaluate`.
- Every dummy argument has an explicit `intent`. Most existing dummy
  arguments comply; new and changed procedures have no exceptions.
- Declarations come at the start of the scope.
- No `save`d state or module variables for data that belongs to one particle;
  per-particle state goes into a derived type that the caller owns
  **(new code)**.

### Kinds and numbers

- Real numbers are `real(dp)` with
  `use, intrinsic :: iso_fortran_env, only: dp => real64` **(new code)**.
  Some older modules use `dp = kind(1.0d0)` or `double precision`; convert
  them in a separate T1 pull request.
- Literal constants carry the kind: `1.0_dp`, `0.5_dp`.
- Use `twopi` and other shared constants from `util` instead of new local
  copies.
- Do not use the IEEE intrinsic modules (`ieee_arithmetic`,
  `ieee_exceptions`) in library sources; they make gfortran save and restore
  the floating-point environment around procedure calls. Use `is_finite` and
  `quiet_nan` from `finite_check`. Tests may use the IEEE modules as an
  oracle.

### Derived types

- Type names end in `_t` (all current types do), for example `tracer_t`,
  `field_can_t`, `symplectic_integrator_t`.
- Extend abstract types only, such as `magnetic_field_t` or
  `coordinate_scaling_t`. To reuse a concrete type, hold it as a component and
  delegate **(new code)**.
- Use `allocatable` components and variables instead of `pointer`. Transfer
  ownership with `move_alloc`. Use `pointer` only where the target must be
  shared or associated, for example procedure pointers for field selection.

### Control flow and logical expressions

- Fortran does not guarantee short-circuit evaluation of `.and.` and `.or.`.
  A guard that protects an index, an allocation, an association or an
  optional argument must be a separate `if`:

  ```fortran
  if (present(ierr)) then
      if (ierr /= 0) return
  end if
  ```

  not `if (present(ierr) .and. ierr /= 0) return`.
- Stop on unrecoverable input or configuration errors with `error stop` and a
  message that names the input key **(new code)**. Report recoverable numerical failures
  through an `ierr` argument or an exit code.

### Performance

- No allocation or deallocation inside the per-step integration loop. Allocate
  work arrays once per particle or per thread, as `simple_main` does for the
  trajectory buffers.
- Loop with the leftmost array index innermost.
- No I/O inside the tracing loop except where output is the purpose; guard
  debug output with the `SIMPLE_ENABLE_DEBUG_OUTPUT` preprocessor switch.
- Measure before and after any performance change and report the result in the
  pull request. Output must stay unchanged unless the pull request says
  otherwise.

### OpenMP and accelerators

- Particles are traced in parallel, one marker per loop iteration, in
  `!$omp parallel` regions with explicit `private` and `firstprivate` lists.
  Prefer `default(none)` or `default(shared)` with explicit private lists in
  new regions.
- Module variables that hold per-thread working state are declared
  `!$omp threadprivate`. Prefer passing state in derived types for new code.
- Shared counters use `!$omp atomic`; use `!$omp critical` only for output
  and other non-numerical serialisation.
- Results must not depend on the thread count beyond floating-point rounding.
  Seeded runs with `deterministic = .True.` must be reproducible.
- Procedures called from OpenACC kernels carry `!$acc routine seq` and must
  not use procedure pointers or `class()` polymorphism in the hot path
  (see `orbit_rk_core`). The CPU path must stay unchanged when acceleration is
  disabled.

### Formatting

- Four-space indentation and lines of at most 88 characters, as configured in
  `.fprettify`. Some older files use two-space indentation or longer lines;
  keep the style of the surrounding code in small edits and do not reformat
  whole files in a functional change.
- Format new files with `fprettify -c .fprettify <file>`. fprettify does not
  read `.fprettify` by itself (it looks for `.fprettify.rc`); without `-c` it
  indents by three spaces. In existing files, inspect
  `fprettify -c .fprettify --diff <file>` and apply only the hunks that touch
  your change. Do not use
  `fprettify --enable-replacements` together with a line-length option; it can
  split tokens across continuation lines.
- CI checks changed Fortran lines with `fprettify -c .fprettify`; run
  `python tools/check_changed_lines.py` before pushing.
- One space before an inline comment. Use `!>` for documentation comments on
  modules, types and procedures.

### Comments and names

- Comments state purpose, units, coordinate conventions and non-obvious
  numerics. History belongs in commit messages and `CHANGELOG.md`.
- Angles are in radians and are not wrapped during integration. State whether
  a radial variable is `s`, `sqrt(s)` or `rho` wherever it crosses an
  interface.
- Names are lower case with underscores. Physical symbols may follow the
  notation of the referenced paper when the comment cites it.

### Generated code

The physics kernels in `src/` are handwritten. CMake generates the version
module in the build tree from `src/version.f90.in`. Prototypes under
`concept/` generate kernels at build time from symbolic derivations. If
generated kernels are added to `src/`, check in both the generator and its
output, mark the output as generated and do not edit it by hand **(new code)**.

## Python

Python code is the `pysimple` package (`python/pysimple/`), the tools in
`tools/`, the examples and the tests in `test/python/`.

- CI uses Python 3.12. Start modules with a docstring and
  `from __future__ import annotations` **(new code)**.
- Follow PEP 8 with a line length of 88 characters. CI runs `ruff check` and
  `ruff format` on changed lines in `python/` and `tools/`
  (`tools/check_changed_lines.py`).
- Public functions have type hints and a docstring in plain prose that names
  units, shapes and raised exceptions **(new code)**. The API reference is
  built from the package docstrings with Sphinx (`docs/`).
- `pysimple` is a thin layer over `simple_api` and `params`. Physics belongs
  in Fortran so that `simple.x` and `pysimple` give the same result.
- Use `pathlib.Path` for paths and `numpy` arrays for state data. States are
  `(s, theta, phi, v/v0, v_par/v)` with shape `(5,)` or `(5, n)`.
- Tests use pytest. Use `pytest.importorskip` for optional dependencies and
  skip, not fail, when optional external data is missing.
- Tools that convert external data, such as `tools/gvec_to_boozer_chartmap.py`,
  document the input and output conventions and the file schema they write.

## Shell, CMake and build files

- Shell scripts start with `#!/bin/bash` or `#!/usr/bin/env bash` and are
  executable (checked by the pre-commit hooks in `.pre-commit-config.yaml`).
- New tests are registered with `add_test` and carry labels (`unit`,
  `integration`, `system`, `python`, `slow`, `regression`) so that the `make`
  test targets select them.
- `cmake/dependency_versions.json` pins fortnum, fortplot, libneo and the
  golden-record libneo revision. Normal builds use the libneo pin unless
  `LIBNEO_REF` or `LIBNEO_PATH` is supplied. Dependency changes are tier T4.
