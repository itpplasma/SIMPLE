# First-order guiding-centre switch study

The recovered FortSym-generated map and its Newton inverse now build against
SIMPLE's native dependency graph. From the SIMPLE root:

```sh
make test CONFIG=Fast TEST='(study_canonical_switch|test_switch_field_jets_)' VERBOSE=1
python concept/fo_gc_derivations/D_canonical_switch/scripts/run_study.py
```

The full study uses all three NCSX and G3V3L3Fi volumes, with the original
tokamak controls, seven epsilon scales, five gyroperiods and temporal refinement
at the finest two scales. Raw data, inputs, native log and fits go under ignored
`build/study/canonical-switch/`; the historical seven-case results remain in
`results/`. The runner invokes a registered native diagnostic. It no longer
compiles against copied libraries or executes a detached binary.

The original SPECTRE provider populated `d2b` with zero. That made the analytic
spatial Jacobian approximate, even when the inverse converged. An independent
Cartesian field-difference oracle rejects that implementation on both
stellarators. The provider now computes logical magnetic-field second jets,
subtracts the geometry-curvature chain-rule term and transforms both derivative
indices with the inverse basis. This is a research derivative provider, with
fixed difference steps inside one smooth owning volume; no production path
uses it. Native field-Hessian checks and the existing switch diagnostic pass
(3/3, 2.55 s). Maximum scaled Hessian error is 2.62e-6 across the six volume
samples; the independent Cartesian spacing difference is at most 7.90e-6.
The threshold was fixed before the repair at 2e-4.

Round-trip identity and O(epsilon^2) magnetic-moment defects alone do not prove
a canonical switch. Independent two-form pullback, phase/frame/normalization,
the extended action-angle GC model, actual sheet entry/exit and long-run loss
validation remain to be completed. These sources are concept-only.

## Completed local checks, 2026-10-14

The ten-case study is retained in `results/native/`: all three NCSX and G3V3L3Fi
volumes, plus four tokamak controls. Finest-four magnetic-moment defect orders
are about 1.995–2.001; zeroth-order and energy-matching inverses remain first
order. Maximum scaled round-trip errors are 2.90e-13 (GC), 2.04e-13 (particle)
and 3.94e-13 (SIMPLE variables). The absolute FO canonical-chart error is
6.65e-11, but relative error is 6.87e-14; it is not an absolute 1e-12 claim.
The largest magnetic-moment defect change on temporal refinement is 0.0328%.
Runtime sources are checked against `f06dbc7` in the retained provenance.

Independent two-form tests compare the generated Jacobian pullback of the
particle form with dGamma for
`Gamma=(A/eps+u*b-eps*mu*Rstar).dX-eps*mu*dphi`,
`Rstar=R+curl(b)/2`, in the same fixed frame and normalization. Analytic fields
include nonzero magnetic Hessians; both stellarators cover all three owning
volumes. Selector 1 defects have fitted order 1.995–2.002; selector 0 supplies
the first-order negative control. Independent finite differences of the actual
map agree with its analytic Jacobian to scaled error at most 8.00e-9. These
are numerical local scaling checks, without production units/phase wiring.

The optional exact checker and regeneration test use an existing owning
FortSym CMake build, leaving ordinary builds CAS-free:

```sh
make test CONFIG=Fast TEST='(test_switch_symbolic_identities|test_switch_regeneration)' \
  VERBOSE=1 FLAGS='-DFORTSYM_BUILD_DIR=/path/to/fortsym/build'
```

The checker passes 132/132 exact identities with zero probes. This includes
cotangent-lift identities, generic finite Lie-series identities and the local
linear, divergence-free field family with symbolic gradient/curvature/shear
parameters. It does not prove a generic nonlinear SPECTRE field identity.
Regeneration reproduces the committed kernel byte for byte. The physical
one-form conventions follow [Tronko–Brizard 2015](https://arxiv.org/html/1511.01005v1)
and [Brizard–Tronko 2016](https://arxiv.org/html/1606.06534v1); the numerical
pullback oracle implements dGamma independently.

The recovered checker initially omitted the local Lie-generator initialization,
replaced perpendicular speed before taking mu derivatives, and substituted
epsilon=0 into unsimplified differentiated expressions. These were checker
construction defects, not refutations of the physical map. Keeping the actual
mu dependence until the coordinate substitution, normalizing coefficients
before taking the limit and declaring the existing mu>0 chart domain closes
those checks. No physical kernel, bound or identity order was weakened.
Production integration remains open: extended canonical action/phase, gyrogauge
connection tables, matching to GC variables, both sheet switches and long-run
loss controls still need implementation and verification.
