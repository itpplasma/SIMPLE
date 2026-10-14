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
