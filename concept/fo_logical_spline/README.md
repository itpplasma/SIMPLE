# Logical-coordinate splined full orbit: stopped W4 prototype

Recovered from `SIMPLE-w4`, branch `feat/spectre-fo-logical-spline`, based on
`45c0598`, after the Claude session stopped on 2026-10-14. See `PLAN.md` for
the programme context, dependencies and acceptance criteria.

The original untracked `gen/`, `src/` and benchmark source are preserved here.
`pending-integration.patch` preserves every tracked edit, including the proposed
replacement of the Cartesian hybrid path. These changes are not active in
production or registered with CMake: the worker stopped before building and
validating the new path. Apply the patch only in a separate implementation
branch after copying the source directories back to the repository root.

The intended field provider splines R, Z, A_theta and A_zeta separately in
each volume and derives the metric from the embedding. The generated canonical
RHS uses those spline jets. The generator reported 11 exact identities in the
old log. Its script still assumes the original root-level `gen/` and
`src/generated/` layout; it is a preserved draft, not a finished reproduction
recipe for this relocated directory.

No speedup or loss-curve agreement is established. Historical Cartesian
baseline timings were 129–138 microseconds/step near rho=0.95 and
88–92 microseconds/step near rho=1.95 on a heavily loaded shared host.
Axis handling, interface gauge continuity, crossing-and-return detection and
equivalence to the sharp-field reference remain acceptance requirements.
