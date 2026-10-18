# Contributing to SIMPLE

Contributions are welcome as bug reports, questions, documentation fixes and
pull requests. This file describes how to build and test SIMPLE locally, how
changes are reviewed and what a pull request should contain. SIMPLE is
maintained by a research group as part of its scientific work. There is no
guaranteed response time for issues or pull requests.

By contributing you agree that your contribution is distributed under the
[MIT License](LICENSE) of this repository.

## Reporting bugs and asking questions

Use the [GitHub issue tracker](https://github.com/itpplasma/SIMPLE/issues) and
choose the matching template:

- **Bug report** for crashes, wrong results, build failures or regressions.
- **Feature request** for new capabilities or changes of behaviour.
- **Question** for usage, input settings and interpretation of output.

A useful bug report contains:

- the SIMPLE commit (`git rev-parse --short HEAD`) or release tag,
- the build command and options (`CONFIG`, `FLAGS`, `LIBNEO_REF`, CUDA or
  OpenACC settings),
- compiler and operating system (`gfortran --version`, `uname -a`),
- the `simple.in` and a description or link of the field input (VMEC `wout`
  file, chartmap, coils),
- expected and actual behaviour, with the relevant part of the output,
- the smallest reproducer you could find: few particles, short `trace_time`,
  `deterministic = .True.`.

Search existing issues first. Do not attach unpublished or confidential
equilibria; describe how a similar case can be produced instead.

## Building

SIMPLE is built with CMake and Ninja through the top-level `Makefile`. Run all
commands from the repository root. The requirements and system packages are
listed in [README.md](README.md#building).

```bash
make                 # Release build in build/
make CONFIG=Fast     # preferred profile for development and physics runs
make CONFIG=Debug    # debugging; also Profile
make clean           # remove build/
```

When changing the profile of an existing build, run `make reconfigure
CONFIG=Fast` (or `CONFIG=Debug`) before `make`. The test targets below, except
`test-all`, reconfigure automatically.

The build fetches libneo, fortplot and fortnum through CMake `FetchContent`.
Select a libneo revision with `make LIBNEO_REF=<branch|tag|sha>` or use a local
checkout with `make LIBNEO_PATH=/path/to/libneo`. The `fpm.toml` is
experimental and is not used for development or testing.

Several tests in every test target, including `make test-fast`, are written in
Python and need `numpy`, `netCDF4` and the other packages in
`requirements.txt`. Create the repository-local environment once and activate
it before running tests:

```bash
./setup-venv.sh
source .venv/bin/activate
```

The default setup also builds and installs `pysimple` and selects OpenBLAS.
Install OpenBLAS first. With Homebrew, make its installation visible to CMake
and select GNU C/C++ as well as GNU Fortran for OpenMP. Replace `16` below
with the installed GCC major version. Select compilers before the first
configuration, or run `make clean` when switching toolchains.

```bash
brew install openblas
CMAKE_PREFIX_PATH="$(brew --prefix openblas)" \
    CC=gcc-16 CXX=g++-16 FC=gfortran ./setup-venv.sh
```

Use `./setup-venv.sh --no-pysimple` to install only the Python dependencies.

Without this environment the Python tests fail with `ModuleNotFoundError`; use
`make test-nopy` to exclude tests labelled `python`. Some remaining test
wrappers still use the Python standard library.

## Running tests

Tests are registered with CTest and run through `make`. Do not call test
executables or `ctest` directly. Except for `test-all`, the targets below
rebuild with deterministic floating point (`SIMPLE_DETERMINISTIC_FP=ON`) first.
`test-all` uses the current build configuration.

| Command | Scope |
| --- | --- |
| `make test-fast` | Unit and integration tests without the `slow`, `regression`, `performance` and `scalability` labels. Run before every push. |
| `make test` | All tests except `regression`, including slow and Python tests. |
| `make test-nopy` | As `make test`, excluding tests labelled `python`. |
| `make test-regression` | Tests labelled `regression`, including golden-record comparisons. |
| `make test-golden` | Golden-record tests only (label `golden_record`). |
| `make test-all` | Every registered test. |

Select tests by name with `TEST=<regex>`, for example
`make test-fast TEST=test_finite_check`. The value passes through `make`, so
write a regex anchor `$` as `$$`. Set `VERBOSE=1` for verbose test output;
otherwise CTest shows output on failure. Add `CONFIG=Fast` to any target to
use the Fast profile.

The golden-record tests clone and build `main` in the build tree, run the cases
in `test/golden_record/*/simple.in` with both builds and compare the outputs.
The first run therefore needs network access and several minutes.

## Tests to add

Every change of behaviour comes with a test that fails without the change.

- **Unit tests** (`test/tests/test_*.f90`, label `unit`) check one procedure or
  module against an independent oracle: an analytic result, a known identity,
  a roundtrip (`ref -> integ -> ref`), a symmetry, or a convergence rate.
  Use them for numerical helpers, coordinate transforms, interpolation and
  input parsing.
- **Integration and system tests** (labels `integration`, `system`) run
  several components or `simple.x` on a small case and check physical
  invariants, for example bounded energy and momentum or agreement between
  field representations.
- **Python tests** (`test/python/test_*.py`, pytest, label `python`) cover
  `pysimple`, the tools in `tools/` and the agreement between `pysimple` and
  `simple.x`.
- **Golden-record cases** (`test/golden_record/<case>/simple.in`) protect
  end-to-end output against unintended change. Add a case when a new field
  type, integrator or output path becomes part of the supported workflow.

Register new tests in `test/tests/CMakeLists.txt` or
`test/python/CMakeLists.txt` with the appropriate labels. Mark tests that run
longer than about a minute as `slow`. Tolerances must follow from a physical or
numerical scale, not from the observed error. A test must not depend on data
that is only available on one machine; missing optional data should skip the
test, not fail it.

A change of golden-record output is a scientific change. The pull request must
state the affected cases, the expected sign and magnitude of the change and the
reason, as described in [docs/review-policy.md](docs/review-policy.md#golden-records).

## Pull requests

1. Create a topic branch from current `main`. Keep one purpose per branch.
2. Make small commits; stage files explicitly.
3. Run `make test-fast` locally, plus the larger targets that your change
   touches (`make test` for Python changes, `make test-regression` for any
   change that can move numerical output).
4. Open a pull request against `main` and fill in the
   [pull request template](.github/pull_request_template.md). Open it as a
   draft while work is in progress; CI runs on pull requests that are ready
   for review.
5. Address review comments with new commits. Pull requests are merged by a
   maintainer with a merge commit.

The full rules are in [docs/review-policy.md](docs/review-policy.md). In
summary:

- Each pull request has one risk tier, from T0 (documentation) to T4
  (parallelism, GPU, dependencies, CI), chosen by the highest-risk part of the
  change. The tier sets the depth of review and the evidence required.
- Each pull request carries exactly one `tier/T0` ... `tier/T4` label and one
  size label (`size/XS` <= 10 changed lines, `size/S` <= 100, `size/M` <= 300,
  `size/L` <= 1000, `size/XL` above). Maintainers add the labels if you cannot.
- The template states the correctness contract: intended change, behaviour
  that must not change, coordinate and unit conventions, numerical invariants,
  tests, golden-record impact and failure modes. A testable claim is evidence;
  a summary is not.
- Split mixed work by risk: feature logic, caller integration, plotting,
  test-performance changes and build or dependency changes go into separate
  pull requests.

For pull requests ready for review, CI builds SIMPLE with deterministic
floating point on Ubuntu and runs the fast, slow and regression tests, a CGAL
build and the documentation build. Follow the
[coding standard](docs/development/coding-standard.md); the current workflows
do not run a formatting or lint check. The workflows in `.github/workflows/`
are authoritative for the current checks.

## Commit messages

Recent history uses a lower-case type prefix and an imperative summary of at
most about 72 characters:

```text
fix: drop the RK momentum gate on accepted symplectic steps
perf: drop the IEEE intrinsic modules from SIMPLE sources
docs: record libneo #427 spline paths and the Boozer GC speedup
```

Common types are `feat`, `fix`, `perf`, `refactor`, `test`, `docs`, `build`,
`ci` and `research` (measurements and studies under `concept/`). Use the body
to explain why the change is needed and what evidence supports it. Reference
issues as `#123`.

## Documentation

Update documentation in the same pull request as the code:

- Input keys: every `&config` key with its default and meaning in
  `examples/simple_full.in` (checked by `test_config_reference`); guidance on
  choosing values in `DOC/config.md`.
- Coordinates, fields, `magfie` modes, splines and canonical coordinates:
  `DOC/coordinates-and-fields.md` must describe the current state after any
  such change.
- Output files and their columns: `README.md` and the relevant `DOC/` page.
- Python API: docstrings in `python/pysimple/` (rendered by Sphinx from
  `docs/`) and `python/README.rst`.
- Changes that can alter results for the same input, and other notable
  changes: one row in `CHANGELOG.md`, with `Y` in the first column when output
  may differ.

Comments in the source explain intent, conventions and non-obvious numerics,
not the history of the change. All angles in SIMPLE are in radians and are not
wrapped during integration; state units and coordinate conventions explicitly
in new interfaces.

## Coding standard

Fortran and Python conventions are in
[docs/development/coding-standard.md](docs/development/coding-standard.md).
