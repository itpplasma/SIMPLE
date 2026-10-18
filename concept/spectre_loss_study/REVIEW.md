# Requested Claude Sonnet reviews, 2026-10-15

Three bounded read-only Claude Code one-shots used `--model sonnet`: FO
provider/events, GC extrapolation, and research/build integration. Each was
limited to five minutes and twelve tool calls, with no edits, tests or agents.
They reviewed integrated tree `6918892` and the focused production commits
`47d3fb0`, `6970eff` and `914544f`. The CLI model alias was requested; a concrete
backend model version was not captured. These are advisory AI reviews, not
independent human domain/senior approval. No new full-suite pass is claimed.

## FO provider and particle-wall events

Sonnet found no established consequential defect in the audited sources.
It checked units, owning gauges, periodic angles, the regular axis and C4
blend, trial/accepted state, loss versus locate failure, single-volume GC
semantics and shared read-only use of the spline tables.

It raised unverified concerns about far-exterior initial locate roots, tiny
logical/Cartesian discrepancies at a graze, and the adjacent volume's radial
basis at an upper sheet. No failing trigger was supplied. No production change
was made on that basis. In particular, a failed Newton solve must not become
a physical loss without independent evidence that the particle is outside.
Marker 129, GC2 and whole-hybrid accuracy were outside this review's scope.

## GC extrapolation

Sonnet independently checked the three corrected terms: midpoint momentum
uses the endpoint radius rather than the independent midpoint radius; the
Gauss mixed Hessian uses packed entry 3; v_parallel includes its p_phi
increment. It found the closed-Hamiltonian oracle appropriately independent,
with limits from its simple h_phi and A_phi and absence of an axis crossing.

It also identified an adjacent, pre-existing axis-chart concern in
`attempt_midpoint_step` and the multi-volume midpoint loop: a negative endpoint
radius is flipped and theta incremented by pi before extrapolating p_theta
from jets evaluated in the old chart. If a converged step takes that path,
the post-flip displacement is not a local Taylor increment in the old chart.
The review called this established by code reading but ran no crossing case.

Controller disposition: retain as an open implementation concern, not a
confirmed explanation of the measured losses. A smooth polar-field probe
did not produce the required accepted negative-radius endpoint and failed its
own chart-crossing precondition. It is not a RED defect witness and was
removed. Existing axis coverage exercises Euler and does not settle this
midpoint path. Before a repair, retain an accepted midpoint crossing with
fresh-field momentum as an independent oracle, then check both midpoint
call paths. Re-evaluation in the new chart, as Euler already does, avoids
assuming canonical-momentum invariance across arbitrary chart changes.
This concern is outside the three-term PR #535 correction.

## Research/build integration

Sonnet found test-only solvers isolated from production, optional CAS
generation guarded, known RED GC energy tests unregistered, and superseded
W4 cleanup consistent with the requested scope.

Its claimed default-suite regression in `test_fo_symplectic_long` and
`test_fo_symplectic_resolved_long` is unsupported for the current source:

- The cited `fo_symplectic_proof/results.json` records predecessor source
  `89babb30`, SHA256
  `2d49f65d6dbef015de5efaeba75e130fe5b9b6715a45ba50fe36f79549185442`.
  That hash matches the retained `reference_gauss2.f90`.
- Current `spectre_fo_gauss2.f90` hashes to
  `2257b247f459ce61f5ee0c5d9b131ee8541bfae1ad3cf00d4316dd575f141818`.
  The registered energy runs use this implementation and its entry step cap.
- The [integrated before/after evidence](../fo_symplectic_proof/README.md#integrated-beforeafter-evidence)
  records the later million-step controls passing unchanged strict gates,
  including the two predecessor Newton stops, at residual tolerance 1e-12.

Historical failed results remain intact. Test registration, tolerances and
labels were not relaxed in response to this finding. No fresh long run was
performed during review; the disposition relies on source correspondence and
the retained later evidence, not a new claim that every default test passes.
The PLAN wording now explicitly distinguishes known RED GC energy acceptance
from these FO controls. Other integration suspicions supplied no failing
trigger and did not justify changes.
