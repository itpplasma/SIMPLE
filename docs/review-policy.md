# Review Policy

SIMPLE changes are reviewed by contract and evidence. A reviewer must be able
to see what behavior changes, what behavior stays fixed, and which independent
checks support that claim.

## Risk Tiers

Each pull request gets one primary tier. Use the highest tier touched by the
change.

| Tier | Scope | Review | Required evidence |
| --- | --- | --- | --- |
| T0 | Documentation, comments, small build metadata | Light review | Relevant build or docs check |
| T1 | Pure refactor with no intended behavior change | Interface review | Golden records unchanged or bitwise behavior unchanged |
| T2 | Local numerical logic, interpolation helpers, coordinate utilities | Focused expert review | Unit tests plus property or roundtrip tests |
| T3 | Physics, output semantics, wall losses, collisions, coordinate conventions | Domain review | Independent comparison plus golden-record explanation |
| T4 | Parallelism, GPU offload, dependencies, CI, security-sensitive build logic | Senior review | Multi-platform evidence, determinism or statistical checks, performance and security review |

AI-generated code is not its own tier. The tier follows the blast radius.

## Human Review Scope

Reviewers should spend time on the contract:

- Is the problem statement correct?
- Is the coordinate convention explicit?
- Which invariant must remain unchanged?
- Which output may change?
- Is the test independent from the implementation?
- Is the pull request small enough to review?
- Is each numerical tolerance tied to a physical or numerical scale?
- Can old input files or post-processing scripts break?

Formatting, syntax, boilerplate, and generated scaffolding belong to tooling.

## Pull Request Contract

Every nontrivial pull request must fill the project pull request template. The
template asks for:

- intended behavior change,
- behavior that must not change,
- coordinate and unit conventions,
- numerical invariants,
- tests added,
- golden-record impact,
- failure modes considered,
- manual validation.

A summary is not evidence. A testable claim is evidence.

## Golden Records

Golden-record changes are scientific changes unless proven otherwise. A pull
request that changes golden records must state:

- changed cases or files,
- expected sign and magnitude,
- physical or numerical reason,
- whether the old behavior was a bug or the new behavior is a feature,
- reviewer acceptance of the scientific movement.

Do not regenerate golden records without explaining the movement. For stochastic
cases, prefer distribution-level checks over particle-by-particle equality.

## Subsystem Invariants

Coordinate transforms:

- `ref -> integ -> ref` roundtrips,
- periodicity,
- explicit `s`, `sqrt(s)`, and `rho` conventions.

VMEC, Boozer, and chartmap equivalence:

- same physical starting surface,
- same short-time orbit within tolerance,
- same loss fraction statistically for comparable runs.

Symplectic integrators:

- bounded invariants,
- no drift beyond the reference envelope.

Collisions:

- Maxwellian fixed point,
- analytic relaxation of pitch and energy moments,
- seeded stochastic reproducibility.

Wall losses:

- same wall hit in Cartesian and coordinate space,
- no missing footprint cluster for reference cases.

NetCDF output:

- semantic equality of `xend_cart`, `wall_hit_cart`, and `zend`,
- explicit units and coordinate frame.

GPU and OpenMP:

- CPU path unchanged when acceleration is disabled,
- short-run CPU/GPU agreement,
- long-run statistical agreement,
- no CPU performance regression.

## Branch Hygiene

Keep pull requests small enough for one invariant owner to review. Split mixed
work by risk:

- feature logic,
- optimizer or caller integration,
- plot or artifact generation,
- test-performance changes,
- build or dependency changes.

T3 and T4 pull requests should carry adversarial review artifacts when possible:
a missing-test pass, an invariant-critic pass, or a split/reduction pass. These
artifacts inform the reviewer. They do not approve the code.
