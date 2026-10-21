# Contributing to SIMPLE

Bug reports, questions, documentation fixes and pull requests are welcome.
SIMPLE is maintained by a research group; there is no guaranteed response
time. Contributions are distributed under the
[MIT License](https://github.com/itpplasma/SIMPLE/blob/main/LICENSE).

## Issues

Use the [issue tracker](https://github.com/itpplasma/SIMPLE/issues) with the
template for a bug report, feature request or question. Search existing
issues first. A bug report contains:

- the commit (`git rev-parse --short HEAD`) or release tag,
- the build command and options (`CONFIG`, `FLAGS`, `LIBNEO_REF`, CUDA),
- compiler and operating system (`gfortran --version`, `uname -a`),
- the `simple.in` and how to obtain the field input,
- expected and actual behaviour with the relevant output,
- the smallest reproducer: few particles, short `trace_time`,
  `deterministic = .true.`.

Do not attach unpublished or confidential equilibria; describe how a similar
case can be produced.

## Build and test

Dependencies and CMake options are in the
[installation guide](https://itpplasma.github.io/SIMPLE/latest/getting-started/installation.html).
For development, from the repository root:

```bash
make CONFIG=Fast              # build in build/
./setup-venv.sh               # once: .venv with test dependencies and pysimple
source .venv/bin/activate
make test-fast                # before every push
```

Select dependencies on the first configuration with
`make LIBNEO_REF=<branch|tag|sha>` or `make LIBNEO_PATH=/path/to/libneo`;
run `make clean` before changing them or the compiler. To change the profile
of an existing build, run `make reconfigure CONFIG=Debug` before `make`.

Run `make test` as well for Python changes and `make test-regression` for any
change that can move numerical output. The
[testing guide](https://itpplasma.github.io/SIMPLE/latest/development/testing.html)
lists all targets and labels, the oracle and tolerance rules and the
invariants a change must keep.

## Pull requests

1. Branch from current `main`, one purpose per branch. Stage files explicitly.
2. Add a test that fails without the change (testing guide).
3. Follow the
   [coding standard](https://github.com/itpplasma/SIMPLE/blob/main/docs/development/coding-standard.md).
   CI checks formatting of changed lines; run
   `python tools/check_changed_lines.py` (needs `fprettify` and `ruff`).
4. Open the pull request against `main` and fill in the
   [template](https://github.com/itpplasma/SIMPLE/blob/main/.github/pull_request_template.md).
   Open it as a draft while work is in progress; CI runs on pull requests
   that are ready for review.
5. Answer review with new commits. A maintainer merges.

Before marking a PR ready, compare the final diff with the issue and PR
contract. Remove unrelated work, redundant tests/docs, transient generated data
and accidental large result sets. Squash fixups or reorganize commits into a few
meaningful steps; preserve rejected scientific approaches only when their
history is useful evidence. Do not rewrite shared, stacked or queued branches
without coordinating with their users.

Each pull request carries one risk tier, set by its highest-risk part. The
tier sets the depth of review and the evidence required:

| Label | Scope | Evidence |
|---|---|---|
| `tier/T0` | documentation, comments, small build metadata | docs or build check |
| `tier/T1` | refactor without intended behaviour change | golden records or bitwise output unchanged |
| `tier/T2` | local numerical logic, interpolation, coordinate utilities | unit tests plus property or roundtrip tests |
| `tier/T3` | physics, output meaning, wall losses, collisions, coordinate conventions | independent comparison and golden-record explanation |
| `tier/T4` | parallelism, GPU, dependencies, CI, security-sensitive build logic | several platforms, determinism or statistical checks, performance |

It also carries one size label for the changed lines: `size/XS` (≤ 10),
`size/S` (≤ 100), `size/M` (≤ 300), `size/L` (≤ 1000) or `size/XL`.
Maintainers add labels if you cannot.

The template asks for the intended change, behaviour that must not change,
coordinate and unit conventions, numerical invariants, tests, golden-record
impact and failure modes. State testable claims; a summary is not evidence.
Split mixed work by risk: feature logic, caller integration, plotting,
test-performance changes and build or dependency changes go into separate
pull requests.

A change of golden-record output is a scientific change. The pull request
names the affected cases, the expected sign and magnitude, the physical or
numerical reason and whether the old behaviour was a bug. Do not regenerate
reference output without that explanation.

## Merging

- Pull requests, including every layer of a stacked PR, run the fast CI tier.
  Label a PR `full-ci` to run the full tier on it before queueing.
- Review a stack layer by layer. When it is ready, collapse it into one PR into
  `main` (retarget the top layer to `main`); the lower PRs then show as merged.
- Click "Merge when ready" (or `gh pr merge --merge --auto`). The PR enters the
  merge queue, the full CI tier runs on the exact merge result, and on success
  it lands on `main` as a merge commit. `main` itself is not re-tested on push.
- A branch need not be up to date with `main`. Merge `main` into a PR only
  when it conflicts.

The merge queue uses `SQUASH`. Standalone PRs land as one commit; stacks
land one squashed commit per layer. Record the lower layer's old head before
merging; rebase the upper layer with
`git rebase --onto origin/main <old-lower-head>` and retarget it to main.
Recheck its final diff and affected tests after rebasing.

The required repository-hygiene workflow scans the PR commit range for secrets,
files over 1 MiB, data/archive types and aggregate additions over 5 MiB. The
[pinned shared checker](https://github.com/lazy-fortran/.github/blob/main/REPO_HYGIENE.md)
also runs offline on staged changes through the prompts global hook, before
this repository's pre-commit hooks. Install `pre-commit` and its existing quick
mechanical checks before committing. Reuse valid test evidence; builds,
numerical tests and optional AI review do not belong in ordinary commit hooks.

Artifact exceptions require a tracked `.repo-hygiene-allow`: repository-relative
path glob, TAB, one-line provenance/reason. There is no environment override.
Never commit non-public machine, experiment or collaborator data as fixtures
(ITER/TC24/IMAS, AUG shots, JT-60SA inputs included). Use synthetic or published
data with provenance. Keep campaign results externally with checksums.
Heavy CI runs only in public repositories. Private repositories validate via
cluster jobs, with job IDs, revision, commands and results recorded in the PR.

## Commit messages

Use a lower-case type prefix and an imperative summary of at most about 72
characters, for example `fix: drop the RK momentum gate on accepted symplectic
steps`. Types: `feat`, `fix`, `perf`, `refactor`, `test`, `docs`, `build`,
`ci`, `research`. The body says why the change is needed and what supports it;
reference issues as `#123`.

## Documentation

Update documentation in the same pull request as the code:

| Change | Update |
|---|---|
| `&config` key or default | `examples/simple_full.in` (checked by `test_config_reference`); usage guidance in `DOC/config.md` |
| coordinates, fields, `magfie` modes, splines, canonical coordinates | `DOC/coordinates-and-fields.md` |
| modules, data flow, extension points | `docs/development/architecture.rst` |
| output files or columns | `docs/user-guide/outputs.rst` |
| Python API | docstrings in `python/pysimple/` and `python/README.rst` |
| results may differ for the same input | one row in `CHANGELOG.md`, `Y` in the first column |

Source comments explain intent, units, coordinate conventions and
non-obvious numerics, not history. Angles are in radians and are not wrapped
during integration.
