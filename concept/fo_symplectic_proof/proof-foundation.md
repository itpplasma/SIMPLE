> Historical proof foundation recovered from `a27d1ec`. Paths and statements
> about the then-shipped experiment refer to that revision. The current
> concept-only layout and verification limits are in `README.md` and `PLAN.md`.

# SPECTRE full-orbit Gauss2: proof chain and executable boundary

This is **not yet a gapless proof of the production hybrid tracer**. It supplies
an exact scalar-generated replacement kernel, a Lean proof for a regular exact
stage branch, and a compiled oracle against a frozen copy of the requested FO
source. Field-provider verification, arithmetic bounds, and the hybrid map have
explicit remaining obligations. Failed experiments are retained in the results.

| Link | Status | Evidence and remaining obligation |
|---|---|---|
| L1: scalar Hamiltonian | MACHINE-CHECKED / OPEN | One scalar in `../../gen/gen_fo_hamiltonian.f90`; A/metric provider is the frozen implementation, with finite derivative checks, not a universal proof of its jets. |
| L2: generated derivatives and emitted source | MACHINE-CHECKED | `symbolic.json`: 257/257 exact identities, zero probes; `monitor.json`: 6/6. `test_fo_symplectic_generated` regenerates and compares all six emitted leaves. |
| L3: exact-stage symplecticity | PROVED | `lean/Gauss2.lean`, theorem `Gauss2.gauss2_symplectic`, registered `test_fo_symplectic_lean`. Parameters and h are fixed; differentiability of the exact stage branch is explicit. |
| L4: nonexact and switching links | BOUNDED / OPEN | Conditional defect inequalities below; finite experiments do not certify their constants. State-selected h, sharp events, entry parameter selection, and GC↔FO maps are outside L3. |
| L5: compiled map | MACHINE-CHECKED / OPEN | Registered monitored NCSX and second-stellarator tests; numerical outcomes and solver limits below. Tokamak is only a short sanity smoke check. |

## L1: what Hamiltonian and field data are actually used

Physical canonical coordinates are q=(rho,theta,zeta), with angles in radians
on an unwrapped cover, and canonical covector p. With mass m>0, charge q_e,
a nonsingular chart metric g, and a C² covariant vector potential A,

```
H(q,p) = (p-q_e A(q))ᵀ g(q)⁻¹ (p-q_e A(q)) / (2m).
Omega = J = [ 0 I ; -I 0 ],       f = J grad H.
```

The scalar generator uses **covariant metric jets**, differentiates the
adjugate/determinant inverse symbolically, and derives every Hamiltonian gradient,
Hessian, RHS and RHS Jacobian. It accepts arbitrary local C² A/metric jets; the
physical expression is not replaced by a polynomial model. Its domain excludes
zero determinant. Positive-definite metric and positive mass are the physical
application conditions.

The frozen step uses normalized momentum P=e_i·v+A_i/ro0 in cm and H=|v|²/2.
The adapter therefore supplies mass=charge=1 and already scaled A. Physical
units and this normalization must not be mixed when replacing the caller.

`reference_gauss2.f90`, SHA-256 in `reference.sha256`, is an unedited snapshot
of `/home/ert/code/SIMPLE-fo-gauss2/src/spectre_fo_gauss2.f90` at
commit `89babb30c0c42504cccd52ec8cfd7ab2d78d38e0`, branch
`feat/spectre-fo-gauss2`. The pin covers the monitored implementation, and the
regeneration test checks the snapshot's SHA256. That external worktree is read
only. SIMPLE/libneo are built locally with `make CONFIG=Fast`; libneo is pinned
to `d30174fb9a9c6520301463ef2406bd2cbfbcd55b`, also used by the source worktree.
The compiled SPECTRE basis files were byte-identical at the initial handoff.

The provider reads exactly these `spectre_field%data` members:

* `Mvol`, `Lrad`, `Mpol`, `mn`, `im`, `in` select volumes and modes.
* `vol(k)%Ate`, `%Ato`, `%Aze`, `%Azo` are the even/odd Fourier coefficients
  of the covariant A_theta/A_zeta radial expansions; A_rho=0.
* Phase is `im*theta-in*zeta`. Outside the axis volume, radial s is
  `2*(rho-owner+1)-1`; inside it, Zernike sbar=(s+1)/2=rho.
* `Rbc`, `Rbs`, `Zbc`, `Zbs` describe the interface embedding. The axis radial
  factor is rho^m (rho² for m=0); other volumes interpolate linearly in rho.
  The Cartesian embedding is (R cos zeta,R sin zeta,Z), scaled by `M_TO_CM`.
* `gauss2_prepare` accumulates volume gauge corrections mode by mode, refusing
  a non-exact tangential jump beyond `1e-11*scale`. A is scaled by
  `TESLA_TO_GAUSS*M_TO_CM**2/ro0`.

Within width/2 of an integer interface, the provider blends **A and the chart
embedding**, using the compact logistic switch in `switch_der`. It does not
smooth B. The metric is the Gram matrix of the blended embedding derivatives.
Gauge projection within a finite tolerance changes the field representation;
it is not a proof that the original interface data match exactly. These tests
use width=0.06; compact layers must not overlap, and the embedding must remain
nonsingular. Axis singularities, exterior extrapolation, and folded charts are
excluded from the theorem's domain.

There is **no spline evaluation in this FO path**: A is evaluated from its
Chebyshev/Zernike–Fourier representation. If a later adapter uses splines, it
must evaluate one A spline and differentiate that same spline; independently
interpolated B or derivative tables do not establish this link.

The actual frozen step integrates the autonomous Hamiltonian

```
K = t(q)*(H-E0),
t(q) = product_k [1-(1-gmin)*exp(-((rho-k)/length)^2)]
       * sqrt(m_ref/(m_ref + g^{rho rho}(q))).
m_ref = (length/(16*h_entry*sqrt(2*E0)))^2 > 0.
```

E0, width, ro0, length, gmin, m_ref and h must be held fixed when
differentiating the map. The physical clock advances by h times the Gauss quadrature of t.
The metric reference uses the entry step, even when a prescribed smaller step
is tested; changing it per step changes the Hamiltonian. Conservation of K does
not imply exact conservation of H: H-E0=K/t, so energy error can be amplified
by small t. Symplecticity alone does not exclude large
physical energy excursions or inadequate resolution of a thin layer.

## L2: generator, interface, and exact checks

`../../gen/gen_fo_hamiltonian.f90` builds H once with symbolic A and symmetric
g functions. FortSym differentiates it, replaces the local derivatives by
independent jet inputs, and emits H, grad H, Hess H, K, J grad K and D(J grad K)
together. The inverse radial metric jet, the metric-monitor factor and its
product with the rho monitor are differentiated from scalar expressions.
The replacement stepper assembles these generated jets without handwritten
Hamiltonian or monitor derivatives.

Jet order is `(value,rho,theta,zeta,rhorho,rhotheta,thetatheta,
rhozeta,thetazeta,zetazeta)`. Metric component order is
`(11,12,13,22,23,33)`. `fo_eval` maps these arrays to the generated scalar
interface. `generated_stepper.f90` only assembles these input arrays, evaluates
the emitted kernel, and solves the 12 stage equations with LAPACK.

The exact checks refuse randomized probe evidence. They check simplification
of all 86 scalar roots, the original unsimplified gradients, the canonical RHS
identity, and Hessian symmetry. They also parse the **actual emitted Fortran**,
expand its straight-line temporaries, and compare every one of its 86 outputs
with its symbolic root. The emitted rho factor, jet product, inverse radial
metric and full monitor leaves are also parsed back and exact-checked. This
proves source-expression correspondence over the stated real-arithmetic domain; it is not a proof of compiler correctness or libm
rounding. Independent energy/RHS and RHS/Jacobian finite differences check the
compiled result, including the supplied field jets.

Normal SIMPLE builds compile committed generated source and need no CAS.
Regenerate from the SIMPLE root with `bash gen/regenerate_fo.sh`, optionally
setting `FORTSYM_BUILD` to an existing FortSym CMake build. Physics stays here;
no FortSym source was changed and no missing generic capability was needed.

## L3: the formal theorem and its hypotheses

`Gauss2.butcher` checks b_i a_ij+b_j a_ji-b_i b_j=0 exactly, with
b=(1/2,1/2), a_11=a_22=1/4, a_12=1/4-sqrt(3)/6,
a_21=1/4+sqrt(3)/6. The algebra works with the actual real sqrt(3), without
replacing it by a decimal.

`quadratic_defect` proves the RK defect of the quadratic invariant
Q(u,v)=omega(u,v) on the doubled tangent space. `symplectic_tangent` cancels
that defect using the Butcher identities. No invariant-preservation or
symplecticity axiom is assumed. `hamiltonian_linearization` obtains
omega(DF u,v)+omega(u,DF v)=0 from C² H, the Hamiltonian identity
omega(F,v)=DH[v], and Mathlib's symmetry of the second Fréchet derivative.
Finally, `gauss2_symplectic` differentiates the exact stage equations and proves
that the derivative of the resulting step preserves omega.

The theorem requires an exact **differentiable branch** of stage solutions.
An exact solution at one input does not imply a unique differentiable map.
To apply it to a solver neighborhood, prove that its stage Jacobian is
invertible and select the corresponding local branch; this branch-construction
obligation is not discharged by a successful Newton iteration at one point.

Lean is pinned to 4.16.0 and Mathlib to
`a6276f4c6097675b1cf5ebd49b1146b735f38c02` in the project manifest.
All printed theorem axiom lists are exactly
`propext`, `Classical.choice`, and `Quot.sound`. There is no `sorry`, custom
axiom, or assumed general symplecticity theorem. The registered test checks
`Gauss2.lean` using the pinned toolchain.

## L4: bounds that do and do not follow

Let the computed smooth stage branches satisfy
Y_i=z+h sum_j a_ij F(Y_j)+r_i(z). Differentiating and applying the same
quadratic-defect identity gives, for unit tangent vectors u,v,

```
(DPhiᵀ Omega DPhi-Omega)[u,v]
 = -h sum_i b_i { omega(DF_i DY_i u, Dr_i v)
                 +omega(Dr_i u, DF_i DY_i v) }.

||defect|| <= 2 |h| ||Omega|| sum_i |b_i| ||DF_i|| ||DY_i|| ||Dr_i||.
```

Thus the defect is O(residual) **in a C¹ residual norm**, on a regular branch
with bounded stage sensitivities. A value bound |r|<=tol alone is insufficient:
r(z)=tol*sin(z/tol²) has small amplitude and unbounded derivative as tol shrinks.
Newton iteration-count changes can also break differentiability of the
terminated map. The tangent returned by both steppers solves the exact
linearized stage equations at the computed stages, and can be symplectic even
when it is not the derivative of the terminated algorithm. This is why the
finite-difference oracle is essential.

For a numerical Jacobian M=M_exact+E with ||E||<=eta and a symplectic exact map,

```
||Mᵀ Omega M-Omega|| <= ||Omega|| (2 ||M_exact|| eta+eta²).
```

Coefficient error additionally contributes h² sum_ij |c_ij|
||Omega|| ||DF_i DY_i|| ||DF_j DY_j||, with c the actual Butcher defect.
These are conditional quantitative bounds, not certified numerical constants
for this binary. Certifying eta requires roundoff, libm, linear-solve,
conditioning and finite-difference truncation bounds. `make test CONFIG=Fast`
selects SIMPLE's deterministic FP profile (O3, native architecture,
ffp-contract=off); the ordinary Fast profile also needs separate fast-math
analysis. No sound libm enclosure or global C¹ Newton residual bound is claimed.

A prescribed sequence of different h values gives a composition of symplectic
maps. A state-selected h=h(z) generally does not: differentiating it introduces
an extra rank-one term. Freezing h after choosing it from z does not remove that
term from the derivative of the complete selection-plus-step map. The same
qualification applies to E0 and layer/refinement parameters selected at FO
entry. Extended-phase-space constructions need their own proof.

With a C² nonsingular blended A/chart there is one Hamiltonian through the
owner change and no separate sharp crossing map. For an unsmoothed interface,
matching tangential A by an exact gauge shift permits a canonical momentum
translation; an energy-matching reset alone is insufficient. Event-time
selection and any reflection/transmission law must also be shown canonical.
Neither the actual sharp event handler nor its roundoff is proved here.

The GC↔FO entry/exit procedure changes representation and chooses gyrophase
and energy-matched velocity. No canonical map between the relevant phase spaces
has been established. The honest hybrid claim is symplectic fixed-parameter
FO segments (subject to the nonexact bounds), with switching and entry parameter
selection outside that claim. No assertion of hybrid symplecticity is made.

## L5: reproducible compiled oracle and integration

The evidence is from **3D stellarator geometry**. Axisymmetric tokamak output
is only a 1,000-step sanity smoke and is excluded from acceptance evidence.
The registered tests use:

* `test/tests/data/ncsx_freeboundary.h5`: Nfp=3, Ntor=2, 18 nonzero toroidal
  modes; its existing provenance describes a low-resolution field/transition
  fixture rather than a force-minimized equilibrium.
* `${libneo_SOURCE_DIR}/test/spectre/data/spectre_test.h5`: the complete
  G3V3L3Fi fixture, Nfp=5, Ntor=2, 18 nonzero toroidal modes. Its adjacent
  `PROVENANCE.md` pins SPECTRE to `f51638520d2dbdee9b0b4e4174ace5b2349b7460`.
  The source is SIMPLE's already pinned libneo dependency; no generated or
  fabricated geometry is needed. A matching copy was found under `~/data`.

`fixtures.json` pins the field hashes and records the search. The scan of
`test/`, `../spectre/` and `~/data/` covered 12,181 HDF5 candidates. The
`../spectre/tests/compare/*/reference.h5` regression caches contain only partial
outputs and cannot serve as complete field fixtures.

Run from the SIMPLE root, using the native authorities:

```bash
make CONFIG=Fast
# Passing short map, field-jet, energy-smoke and saved-failure diagnostics:
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 make test CONFIG=Fast \
  TEST='^test_fo_symplectic(_resolved|_second|_second_resolved|_newton_boundary)?$' VERBOSE=1
# Exact regeneration/readback and the pinned Lean theorem:
make test CONFIG=Fast TEST='^test_fo_symplectic_(lean|generated)$' VERBOSE=1
# Long experiments, including retained failures described below:
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 CTEST_PARALLEL_LEVEL=4 \
  make test CONFIG=Fast \
  TEST='^test_fo_symplectic(_resolved|_second|_second_resolved)?_long$' VERBOSE=1
```

The passing results are both velocity variants of each stellarator short test,
Lean, regeneration/diff, the saved Newton-boundary test, and both second-
stellarator long tests: **9 passing registered tests**. Both NCSX long tests
return failure for the retained cases below. They are not disabled, marked
expected-success, or accepted with a weaker residual criterion. Initial attempts
to use an incomplete second-fixture regression cache were discarded; final
second-fixture results use the complete pinned field.

The stepper's absolute Newton tolerances are 1e-7, 1e-10 and 1e-12. Each is
measured for both 1 and 1,000 steps, at rho=.65 and rho=.995, with h=.02,
ro0=1e5, width=.06 and q angles (.7,.3). Baseline Cartesian velocity is
(.001,.002,.001); `resolved` multiplies it by 100. The rho=.65 tests use
length=1, gmin=1; rho=.995 uses length=.12, gmin=.1. E0 and m_ref are set
once at the unperturbed start and held fixed in every finite-difference call.
The frozen reference API exposes neither tolerance nor residual, so it prints
-1 for those fields. Its stopping rule differs from the strict replacement.

The map oracle uses fourth-order central differences at spacings
1e-4*max(1,|z_j|) and 5e-5*max(1,|z_j|), with a separately composed implicit
stage tangent. Endpoint-only finite differences omit discarded tangent solves;
the unperturbed map still composes every stage tangent. The frozen provider's
fixed-tolerance sweep runs once, retaining its strict 1,000-step gate.
Norms are Frobenius in the stated canonical coordinates.
The tangent is the derivative of the exact stage equations at the computed
stages, not automatically the derivative of the terminated Newton algorithm.
Strict thresholds remain 1e-7 for the tangent, 1e-5 for map finite differences,
1e-6 for compiled derivative consistency and 1e-8 for endpoint correspondence.
The original strict gates cover all one-step maps and the 1,000-step map at
1e-12 in the slow ``_full`` cases. The fast tier composes 64 steps at the same
strict thresholds; its energy sample counts, field-jet and derivative checks
are unchanged. The additional loose-tolerance 1,000-step sweep is an audit: exceeding
that same strict criterion is explicitly printed and retained as a finding.
The reference comparison at loose tolerance is reported rather than requiring
an algorithm with a different stopping rule to return the identical endpoint.

### Map defects versus Newton tolerance

Values below are maxima over the two starting radii. The last column is the
Frobenius disagreement between the two finite-difference Jacobians, not an
error enclosure.

| Field / velocity | Newton tol | FD defect, 1 step | FD defect, 1000 steps | Tangent defect, 1000 steps | FD disagreement, 1000 steps |
|---|---:|---:|---:|---:|---:|
| NCSX, baseline | 1e-07 | 2.914e-10 | 1.363e-05 | 9.549e-14 | 2.671e-04 |
| NCSX, baseline | 1e-10 | 2.966e-11 | 3.490e-09 | 7.391e-14 | 1.420e-07 |
| NCSX, baseline | 1e-12 | 2.966e-11 | 3.490e-09 | 7.391e-14 | 1.420e-07 |
| NCSX, 100x | 1e-07 | 6.070e-11 | 2.030e-05 | 4.768e-13 | 2.770e-05 |
| NCSX, 100x | 1e-10 | 6.070e-11 | 1.614e-07 | 2.252e-13 | 1.180e-06 |
| NCSX, 100x | 1e-12 | 6.070e-11 | 1.614e-07 | 2.252e-13 | 1.180e-06 |
| G3V3L3Fi, baseline | 1e-07 | 2.748e-11 | 7.389e-09 | 1.555e-14 | 1.290e-09 |
| G3V3L3Fi, baseline | 1e-10 | 2.748e-11 | 1.534e-09 | 2.343e-14 | 1.186e-09 |
| G3V3L3Fi, baseline | 1e-12 | 2.792e-11 | 1.710e-09 | 2.343e-14 | 1.313e-09 |
| G3V3L3Fi, 100x | 1e-07 | 6.049e-10 | 2.118e-08 | 1.486e-14 | 1.353e-08 |
| G3V3L3Fi, 100x | 1e-10 | 2.375e-10 | 2.118e-08 | 1.081e-14 | 1.353e-08 |
| G3V3L3Fi, 100x | 1e-12 | 2.375e-10 | 2.118e-08 | 1.081e-14 | 1.353e-08 |

The 1e-7 NCSX sweeps exceed the strict finite-difference criterion in one
1,000-step map each. Stopping-branch changes and finite-difference resolution
are visible in the large spacing disagreement; these are not evidence for an
exact symplectic terminated map. At 1e-12 the largest NCSX FD defect is
1.614e-7, with spacing disagreement 1.180e-6. Thus the experiment resolves only
an approximately 1e-6 scale there, not a certified 1.614e-7 defect. The largest
second-stellarator 1e-12 FD defect is 2.118e-8. All strict map gates pass.

Independent physical-H checks compare value/RHS/Hessian to the frozen source.
Additional monitored checks evaluate the scalar K independently and compare its
finite-difference gradient and the RHS derivative to the generated kernel.
They use m_ref=g^{rho rho} at the initial test point and E0=.9H to exercise the
full second-monitor terms away from H=E0. Spatial spacing refinement retains
the same 1e-6 derivative threshold; this is finite evidence about the provider,
not a universal proof of its jets.

### Million-step energy envelopes and fitted slopes

Each completed row has 1,000,000 accepted steps at absolute residual tolerance
1e-12. H is evaluated after every step. The fit uses all per-step samples of
H/E0-1 against the fictitious step index, accumulated by centered online
regression. Slopes and standard errors have units of relative energy error per
step. The standard error is the ordinary least-squares value; orbit samples are
correlated, so it is descriptive and cannot certify statistical absence of drift.

| Field / start | min(H/E0-1) | max(H/E0-1) | linear-fit slope / step | OLS standard error / step |
|---|---:|---:|---:|---:|
| NCSX rho=.65 | -9.054535e-12 | 4.705480e-11 | 2.903802e-17 | 3.291132e-20 |
| NCSX rho=.95 | -1.356248e-12 | 5.097212e-11 | 7.392167e-18 | 2.797282e-20 |
| NCSX rho=.995 | -3.912859e-11 | 3.330713e-11 | -4.572588e-17 | 3.643221e-20 |
| NCSX seed 9 | -4.171227e-10 | 9.305615e-09 | -2.367611e-16 | 1.678954e-18 |
| NCSX seed 20 | -7.540495e-06 | 6.341560e-09 | 8.763618e-13 | 5.375115e-15 |
| NCSX seed 21 | -2.833332e-10 | 6.174738e-09 | -1.796832e-16 | 1.367793e-18 |
| NCSX 100x rho=.95 | -4.720567e-08 | 3.774922e-08 | 1.288911e-17 | 1.586817e-18 |
| NCSX 100x rho=.995 | -3.304679e-11 | 2.609513e-11 | 1.309093e-19 | 2.309975e-21 |
| G3V3L3Fi rho=.65 | -2.203582e-11 | 1.737810e-11 | -3.064189e-17 | 2.422828e-20 |
| G3V3L3Fi rho=.95 | -3.389400e-12 | 3.847966e-11 | 3.433857e-17 | 1.421273e-20 |
| G3V3L3Fi rho=.995 | -3.095302e-12 | 2.385225e-11 | 1.110245e-17 | 1.597889e-20 |
| G3V3L3Fi 100x rho=.65 | -1.580958e-13 | 1.061373e-13 | -1.004775e-19 | 1.241903e-22 |
| G3V3L3Fi 100x rho=.95 | -1.116884e-13 | 1.325606e-13 | 3.144996e-20 | 1.816719e-22 |
| G3V3L3Fi 100x rho=.995 | -1.216804e-13 | 7.640555e-13 | -2.358362e-20 | 7.438778e-22 |

The interior and two edge starts on each stellarator use the parameters above.
NCSX seeds 9, 20 and 21 are identifiable from the predecessor validation set:
`g_runs/fixed_{9,20,21}/start.dat` under the original job directory. Their GC
coordinates and pitch are retained in `validation_initial`. The test constructs
a gyrophase-zero FO entry using the production velocity, Larmor displacement,
encoding and frozen-parameter formulas at those seeds (3.5 keV alpha ro0,
speed sqrt(2)). They are not later checkpoints of the hybrid trajectories and
do not reproduce all GC history. These physical-speed seeds already exceed the
(.1,.2,.1) larger velocity used for the synthetic 100x starts. Their full initial
canonical states and parameters are retained in `results.json`.

These runs satisfy the unchanged 1e-5 finite energy-envelope criterion. The
observed slopes are not all zero: seed 20 has a fitted 8.764e-13 per-step trend,
and the very small baseline slopes include accumulated numerical drift.
Bounded observed envelopes do **not** prove an all-time bound, and “no secular
drift” as a universal or statistically certified claim remains OPEN. In
particular the phase-dependent seed-20 envelope spans several parts in 1e6;
the table retains its nonzero fitted slope alongside the extrema.

### Retained failures and their diagnosis

1. **The predecessor's exact larger-velocity start**, rho=.65, h=.02,
   width=.06, ro0=1e5, length=gmin=1, now with the pinned metric monitor,
   refuses step **379060** at rho=**.9888925084970212**, residual
   **.0250120161** against 1e-12. The old run refused step 379044; that
   failure persists with the monitor.
   The frozen reference also fails this state. Independent K/gradient and
   RHS/Jacobian checks are below 1.8e-9 at the final refined spacing.
   A prescribed h/2 step, with every Hamiltonian parameter including m_ref
   unchanged, converges with residual **2.4446e-14**; h/4 gives 4.7413e-14.
   This indicates a nonlinear solve / stage-resolution limitation, without
   evidence of a kernel derivative mismatch. Existence and basin of a regular
   full-h stage branch remain unproved. This setup has gmin=1 and m_ref=162.7604;
   it is not the much more strongly refined production layer configuration.
   Its missing million-step envelope remains OPEN; larger-velocity edge starts
   do complete the requested million steps.
2. **Additional captured marker-20 FO state**, from the pinned source's
   radial-refinement regression, starts at rho=1.0025024863675955 with width
   .0018214223148575103, h=.0375521103275432, length=.05828551407544033,
   gmin=.001953125 and production ro0. The strict stepper refuses step
   **151233**, rho=.999652142971619, residual **1.26511301e-12**.
   The reference accepts it under its different scaled/update stopping rule;
   a prescribed half step reaches residual 1.8004e-13. At successively refined
   finite-difference spacings, the RHS/Jacobian discrepancy falls from
   1.1941e-3 to 4.7025e-6 to **2.5959e-8**. No derivative inconsistency is
   demonstrated. A strict-residual/roundoff-floor limitation is plausible,
   but not certified by these probes. No tolerance or iteration limit was
   changed. This extra million-step experiment is incomplete; the first
   100,000-step block stayed within [-1.745e-8,1.717e-8].

`results.json` retains all final test outcomes, map measurements, completed
energy fits, starting states and failure diagnostics. Recreate it from captured
verbose make-test logs with `python concept/fo_symplectic_proof/collect_results.py
LOG...` (newer reruns replace earlier results for the same registered test).
Finite fixture tests do not establish every marker, chart, smoothing width,
Newton branch or production entry setting.

### Plugging the generated kernel into spectre_fo_gauss2

1. Add `src/generated/fo_kernel.f90` and its include directory to the SIMPLE
   target. Normal builds use committed generated source and need no CAS.
2. Reuse the packing in `generated_eval` to form aj(3,10) and gj(6,10) from
   the **same** spectral A and smooth embedding provider used by H.
3. Aggregate rho factors with `fo_monitor_factor` and `fo_jet_product`.
   Compute `rj(10)` with `fo_radial(gj,rj)`; compute the full monitor jet
   `tj(10)` with `fo_metric_monitor(rj,rhoj,m_ref,tj)`. Hold m_ref fixed at
   FO entry using the pinned caller's formula. The zero/omitted reference
   compatibility case uses only the rho jet, packed into slots 1,2,5.
4. Call `fo_eval(P,aj,gj,1,1,E0,tj,H,K,gradH,hessH,f,df)` and use f and df
   directly in the stage equations. gradH/hessH are physical H derivatives;
   f/df already contain every derivative of K. Remove the handwritten
   transformed-Hamiltonian and inverse-metric derivative assembly.
5. Advance the physical clock using h/2 times the sum of the two stage
   monitor values tj(1). Preserve failure propagation and every frozen entry
   parameter. Selecting h or m_ref from the evolving state changes the map.
6. Close the provider, regular-branch, arithmetic and solver limits before
   claiming gapless binary symplecticity; audit GC/FO selection and switching
   separately. No production dispatch or external worktree was changed.
