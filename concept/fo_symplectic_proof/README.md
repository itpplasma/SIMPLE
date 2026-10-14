# Experimental Gauss2 FO integration

This implementation ships the scalar-generated kernel from proof commit
`a27d1ec50d32e853b3ea0806e4662ed8d1be5b0c`. The proof repository at
`/home/ert/code/SIMPLE-sproof`, branch `concept/fo-symplectic-proof`, is read-only
integration evidence. Its README describes the exact L1–L5 boundary and the
Lean theorem `Gauss2.gauss2_symplectic`.

The six `src/generated/` leaves, scalar generator, and regeneration/diff test
are unchanged from that commit. SIMPLE compiles the committed kernel without
requiring a CAS. The registered regeneration test re-emits and compares all six
leaves and runs the exact symbolic/readback checks (257/257 Hamiltonian and 6/6
monitor identities). Thus the shipped expressions remain tied to L1/L2.

## Map and scope

The canonical variables are q=(rho,theta,zeta) on the unwrapped angular cover,
P=e_i·v+A_i/ro0_bar, with unit normalized mass and charge. The spectral provider
blends covariant A and the chart embedding, forms covariant metric jets, and
packs the same ten-slot jet order used by the generator. `fo_radial`,
`fo_metric_monitor` and `fo_eval` supply the complete derivatives of

```
H = (P-A/ro0_bar)^T g^{-1}(P-A/ro0_bar)/2,
K = t(q)*(H-E0).
```

The physical clock uses h/2 times the sum of the two stage monitor values.
Newton solves for stage increments rather than subtracting two absolute
momenta. Acceptance requires absolute maximum stage residual <=1e-12;
there is no update-size or scaled-residual escape and the 20-iteration limit
is unchanged. Failure leaves the last accepted canonical state intact.

Production selects `spectre_fo_experimental_gauss2=.true.`; RK45 remains the default.
At FO entry, h_entry=0.1*ro0_bar/|B_entry| and radial_speed is the local metric
bound on speed in rho. The layer width is 0.05*ro0_bar/|B_entry|*radial_speed,
length=max(16*h_entry*radial_speed,16*width), and
gmin=min(0.1,width/(256*h_entry*radial_speed)). E0 and
m_ref=(length/(16*h_entry*sqrt(2*E0)))^2 are frozen there. The chosen step is
h=h_entry/2*min(1,width/(gmin*length)), fixed throughout the excursion. An explicit convergence-study
factor can further reduce h without changing m_ref or K.

This prescribed refinement uses only entry data and is decided before the
first stage solve. There is no reject/retry or Newton-driven step selection.
For any fixed set of entry parameters every accepted step is the same fixed-h
Gauss2 map, so a regular differentiable exact-stage branch meets L3. The
complete map that selects E0, h and other parameters from its input is outside
L3: freezing selected parameters does not remove their derivatives from that
complete map. GC↔FO switching is not established canonical. Output-clock
overshoot is carried forward rather than shortening a Gauss step.

L4 remains open for global C¹ terminated-Newton residual bounds, arithmetic
and libm bounds, provider-jet certification and regular stage-branch selection.
Nonsingular C² blended charts, nonoverlapping layers and exact compatible
interface gauges are hypotheses. Axis singularities and exterior extrapolation
are outside the theorem. Observed energy envelopes do not prove all-time
boundedness or a universal absence of secular drift.

## Two NCSX failures and marker-20 diagnosis

The predecessor larger-velocity run refused step 379060 at rho=.9888925,
residual .025012 against 1e-12. The saved-state regression checks the original
full-h solve and prescribed half/quarter steps, without altering any kernel or
Hamiltonian parameter. The integrated experiments freeze the layer-resolution step cap above.
The unrefined synthetic case has length=1, width=.06, gmin=1: its original
metric monitor only bounds the continuous on-shell radial displacement by
length/16. The selected h=.0006 instead bounds that scale by width/32. This is
a resolution guide for the underlying flow, not a proof of a finite-step bound.
For the strongly monitored production seeds, width/(gmin*length)>1, so h is
halved without changing any Hamiltonian parameter.

A constant h/2 trial still fails at step 766151 (residual .056689). A constant
h/4 trial completes the same fictitious duration in 4M steps but reaches an
energy excursion 3.52198e-2. These failed trials establish that Newton success
alone does not resolve this case. At the predecessor's saved Newton state,
H=.02947208 against entry E0=.03, already a -1.76% energy defect; a later
h/2 saved state is also off-shell. Continuing such a state cannot restore the
original energy by a symplectic step, and physical H need not stay constant
under K=t(H-E0) off-shell. The saved-state Newton test retains the original
E0 and m_ref; the matched-duration run starts at the original on-shell state.

The seed-20 trend and the captured marker-20 Newton stop are distinct starts.
The seed is a gyrophase-zero FO entry from the retained GC validation seed;
the captured state starts at rho=1.0025024863675955 with a narrower layer.
All diagnosis runs below are pure FO, with no hybrid switch.

Matched fictitious duration for seed 20:

| h/base | Steps | min(H/E0-1) | max(H/E0-1) | slope/accepted step |
|---|---:|---:|---:|---:|
| 1 | 1,000,000 | -7.54055144e-6 | 6.34142716e-9 | 8.76387773e-13 |
| 1/2 | 2,000,000 | -4.71482044e-7 | 3.96092714e-10 | 2.73971362e-14 |
| 1/4 | 4,000,000 | -2.98047860e-8 | 2.42883491e-11 | 8.19682524e-16 |

The envelope reduction factors are 15.99 and 15.82, consistent with fourth-order
resolution error. At the same first 50,000 base steps, tightening Newton's
residual from 1e-12 to 1e-13 changes the negative envelope from -6.07882895e-6
to -6.07882402e-6 (4.93e-12), six orders below the excursion. The longer 1e-13
probe stops at step 74405 with residual 1.01089276e-13. Resolution therefore
dominates the measured seed-20 error; these measurements do not establish a
global absence of branch switching or certify roundoff constants. Slopes remain
nonzero and describe correlated samples; their OLS standard errors do not
certify absence of drift. Per-step slopes at different h have different units
of elapsed time; multiplying by the number of matched-duration steps gives
trend contributions 8.764e-7, 5.479e-8 and 3.279e-9.

For the captured marker-20 state, the original absolute-stage formulation
stopped at step 151233 with residual 1.26511301e-12. Solving stage increments
at the original h moves the stop to 280914, residual 1.01053541e-12; it does
not eliminate the strict residual floor. Fixed half-size steps complete a
million with envelope [-5.07856635e-9,1.02360653e-9], slope -1.80948040e-16
and maximum residual 9.99988287e-13. The finite derivative probes and the
resolution/tolerance comparison support a numerical resolution/residual-floor
cause rather than a demonstrated generated-kernel mismatch. No tolerance or
energy/defect threshold was relaxed.

## Registered evidence and reproduction

Acceptance evidence uses NCSX and the complete G3V3L3Fi fixture only. Their
hashes and nonaxisymmetric mode metadata are retained in `fixtures.json`.
The original proof used libneo d30174f; the integrated build uses c92b336,
whose only intervening change is in the RK ODE driver. SPECTRE source and the
G3V3L3Fi fixture are unchanged. Tokamak is field/tangent/periodicity smoke and
a 1,000-step sanity run only, excluded from acceptance measurements.

The map oracle calls the shipped implementation. The unedited hashed
`reference_gauss2.f90` is compiled with only its module name changed in the
build tree to compare the predecessor independently. The finite-difference
oracle checks the terminated map as well as its implicit-stage tangent. If the
original spacings 1e-4/5e-5 leave a strict-tolerance 1,000-step map unresolved,
it also evaluates 2.5e-5/1.25e-5, retaining the coarse result and applying the
same 1e-5 defect threshold. The refined spacing disagreement must itself be
<=1e-5. This increases resolution without relaxing the defect gate.
Acceptance thresholds remain 1e-7 (tangent), 1e-5 (map FD), 1e-6 (derivative
consistency), 1e-8 (endpoint correspondence), and 1e-5 (long energy envelope).
The loose-tolerance 1,000-step sweeps retain excursions beyond the strict FD
criterion as audit findings, as in the proof evidence.

```bash
make CONFIG=Fast
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 CTEST_PARALLEL_LEVEL=6 \
    make test CONFIG=Fast VERBOSE=1
# Optional matched-duration and tighter-tolerance diagnoses:
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 CTEST_PARALLEL_LEVEL=4 \
    make test-all CONFIG=Fast TEST='test_fo_gauss2_resolution_' VERBOSE=1
# Regeneration can use an existing FortSym build:
FORTSYM_BUILD=../fortsym/build bash gen/regenerate_fo.sh
```

The diagnosis commands include the intentionally retained base-h captured-state
failure; it is not marked expected-success. `collect_results.py LOG...` records
registered maps, per-marker envelopes/slopes and solver diagnostics. All
million-step fits use every accepted sample with centered online regression.
One million accepted steps cover a shorter fictitious duration under the entry
resolution cap. The separate 2M/4M seed-20 rows compare equal duration, and the
33,333,334-step larger-velocity policy run covers the predecessor's 20,000
fictitious time units (last fixed step overshoots by .0004).

## Integrated before/after evidence

All after rows contain 1,000,000 accepted steps at the unchanged absolute
1e-12 residual tolerance. Before rows are the pinned a27d1ec experiments.
Slopes have units of relative physical-energy error per accepted fictitious
step; h differs, so these rows do not alone compare equal duration.

| Field / start | Before [min,max] | After [min,max] | Before slope/step | After slope/step |
|---|---:|---:|---:|---:|
| NCSX rho=.65 | [-9.054535e-12, 4.705480e-11] | [-7.016943e-12, 3.936673e-11] | 2.903802e-17 | 2.839681e-17 |
| NCSX rho=.95 | [-1.356248e-12, 5.097212e-11] | [-2.460787e-11, 2.952993e-11] | 7.392167e-18 | -2.483563e-17 |
| NCSX rho=.995 | [-3.912859e-11, 3.330713e-11] | [-4.922895e-11, 1.735878e-11] | -4.572588e-17 | -3.826867e-17 |
| NCSX seed 9 | [-4.171227e-10, 9.305615e-09] | [-3.688883e-11, 5.823739e-10] | -2.367611e-16 | -1.568731e-17 |
| NCSX seed 20 | [-7.540495e-06, 6.341560e-09] | [-4.714820e-07, 3.960927e-10] | 8.763618e-13 | 3.101830e-14 |
| NCSX seed 21 | [-2.833332e-10, 6.174738e-09] | [-7.163714e-11, 3.859046e-10] | -1.796832e-16 | -3.493254e-17 |
| NCSX captured 20 | stopped 151233 | [-5.078566e-09, 1.023607e-09] | incomplete | -1.809480e-16 |
| NCSX 100x rho=.65 | stopped 379060 | [-3.929079e-13, 2.151612e-13] | incomplete | -4.670751e-19 |
| NCSX 100x rho=.95 | [-4.720567e-08, 3.774922e-08] | [-8.137935e-14, 3.848033e-13] | 1.288911e-17 | 2.858903e-19 |
| NCSX 100x rho=.995 | [-3.304679e-11, 2.609513e-11] | [-2.573497e-13, 2.932321e-12] | 1.309093e-19 | 8.312073e-19 |
| G3V3L3Fi rho=.65 | [-2.203582e-11, 1.737810e-11] | [-1.936040e-11, 1.433920e-11] | -3.064189e-17 | -1.366039e-17 |
| G3V3L3Fi rho=.95 | [-3.389400e-12, 3.847966e-11] | [-8.636314e-12, 2.965317e-11] | 3.433857e-17 | 2.250276e-17 |
| G3V3L3Fi rho=.995 | [-3.095302e-12, 2.385225e-11] | [-2.688749e-11, 8.728795e-12] | 1.110245e-17 | -2.556340e-17 |
| G3V3L3Fi 100x rho=.65 | [-1.580958e-13, 1.061373e-13] | [-3.319567e-13, 7.260859e-14] | -1.004775e-19 | -3.128081e-19 |
| G3V3L3Fi 100x rho=.95 | [-1.116884e-13, 1.325606e-13] | [-2.045031e-13, 3.197442e-14] | 3.144996e-20 | -2.258514e-20 |
| G3V3L3Fi 100x rho=.995 | [-1.216804e-13, 7.640555e-13] | [-3.819167e-14, 1.662892e-12] | -2.358362e-20 | -3.799590e-20 |

The original-start larger-velocity NCSX run also completes **33,333,334**
steps at h=.0006, covering the original million-step fictitious duration.
Its envelope is [-6.931390e-07, 9.077989e-08], slope 5.17708984e-19,
OLS standard error 2.01006712e-20 per accepted step, and maximum
residual 9.99907188e-13. It has no Newton failures.

Strict map maxima below cover the two start radii, at tolerance 1e-12,
including one-step and 1,000-step maps. These are finite observations.

| Field / velocity | Before FD defect | After FD defect | After tangent defect | After FD spacing disagreement |
|---|---:|---:|---:|---:|
| NCSX | 3.489900e-09 | 1.416800e-09 | 2.120300e-14 | 5.411300e-08 |
| NCSX 100x | 1.613600e-07 | 1.589700e-07 | 2.825000e-13 | 5.300300e-06 |
| G3V3L3Fi | 1.709900e-09 | 1.018200e-09 | 1.172600e-14 | 8.045300e-10 |
| G3V3L3Fi 100x | 2.117600e-08 | 1.352100e-08 | 1.217200e-14 | 1.031600e-08 |

The NCSX 100x FD value is below its spacing disagreement and is not a
certified defect enclosure. All unchanged strict gates pass. The selected
entry policy eliminates both retained Newton stops in the requested on-shell
runs; the saved-state check also accepts the new fixed h with the original
E0 and metric reference. Coarse optional diagnosis profiles still retain
their failures.

`make CONFIG=Fast` succeeds. The final full `make test CONFIG=Fast` passes
**96/96**, including generation/diff, in 549.73 seconds. The matched-duration
registered run and its associated short oracle pass 2/2 in 1564.13 seconds.
No ctest command or fpm was invoked directly.

The remaining open items are L4 (provider/domain certification, regular branch
selection, C¹ residual and arithmetic bounds), the noncanonical hybrid switch
and entry/output selection, universal all-time energy bounds, and performance
of strongly refined FO excursions. Gauss2 FO remains experimental.
