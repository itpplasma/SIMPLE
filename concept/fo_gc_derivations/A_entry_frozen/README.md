# A: entry-frozen parameters for the SPECTRE FO Gauss2 map

Stopped early when the scope changed (Gauss2 FO is concept-only). The Lean file
`concept/fo_symplectic_proof/lean/EntryFrozen.lean` now contains **only the
theorems that check** (0 errors, every `#print axioms` exactly
`propext, Classical.choice, Quot.sound`, no `sorry`). The declarations that
failed to elaborate (`selection_hasFDerivAt`, `locally_constant_derivative`,
and the whole `Plane` section: switch example and the H=p^2/2 counterexample)
were moved unchanged to `EntryFrozenDraft.lean.txt` in this directory. That
draft is not part of the Lean project and is **not machine-checked**; statuses
below that mention `sorryAx` refer to it.

Check (about 1 min):
`cd /home/ert/code/SIMPLE-derive && lake +leanprover/lean4:v4.16.0 -d concept/fo_symplectic_proof/lean env lean concept/fo_symplectic_proof/lean/EntryFrozen.lean`

## Claims and status

1. **Fixed-lambda composition: MACHINE-CHECKED.** `symplAt_comp`,
   `compose_symplectic` (by induction over N) and `frozen_composition`. If one
   fixed lambda fixes K, F and h, and each step has an exact differentiable stage
   branch at its orbit point, then the N-step composite is symplectic. This uses
   `Gauss2.gauss2_symplectic`.
2. **Leaf/Poisson statement: MACHINE-CHECKED.** `extend_leaf` and
   `extend_on_leaf`: G(z,lam)=(Phi_lam z, lam) maps each leaf to itself, and its
   restriction to a leaf is Phi_lam. `extend_hasFDerivAt` gives DG.
   `leaf_poisson_iff`: G is a Poisson map of the Casimir-extended bivector
   exactly when D_zPhi preserves the leaf bivector. D_lamPhi does not enter.
   `leaf_poisson_of_symplectic`: if A preserves a nondegenerate omega and is
   surjective, A preserves the sharp-defined bivector. Surjectivity is a
   hypothesis here; in finite dimensions it follows from nondegeneracy.
3. **Selection obstruction: PARTIAL.**
   - MACHINE-CHECKED: `selection_defect` and `selection_defect_symplectic`. These
     are exact identities with C = D_lamPhi∘Dsigma and A = D_zPhi:
     omega(Mu,Mv)-omega(u,v) = omega(Au,Cv)+omega(Cu,Av)+omega(Cu,Cv).
     `locally_constant_selection` is also checked: if sigma is locally constant
     near x, Psi is symplectic at x.
   - NOT CHECKED: `selection_hasFDerivAt`, which gives
     DPsi = D_lamPhi∘Dsigma + D_zPhi. It failed because the `HasFDerivAt.comp`
     argument needs a `show`, which was not added before the stop. This also
     leaves `locally_constant_derivative` unchecked. By hand it is the chain rule.
   - NOT CHECKED: the switch example (`switch_not_continuous`,
     `switch_not_symplectic`, `switch_symplectic_off`). A `ContinuousAt.comp`
     unification failed and a typeclass problem got stuck. The mathematics is
     elementary: the map q+s(q) jumps by 1 at q=0.
4. **Counterexample H=p²/2: PARTIAL.** `free_stage₁` and `free_stage₂` are
   MACHINE-CHECKED: with the actual Gauss2 coefficients, the stage values are
   Y_i = (q + h(a_i1+a_i2)p, p). The rest depends on `sorryAx`: `free_step`,
   which reduces the step to (q+hp, p), failed on one leftover `True ∨ h=0`
   goal, and `free_hasFDerivAt` failed because `HasFDerivAt.pow` does not exist
   in this Mathlib. The theorems `free_ham`, `free_fixed_symplectic`,
   `selected_hasFDerivAt`, `selected_pairing` and `selected_not_symplectic` are
   therefore unchecked. Hand result: Psi(q,p) = (q+h(q)p, p), and
   omega(DPsi e1, DPsi e2) = 1 + h'(q)p. For h = h0+q and p = 1 this is 2, not 1.
5. **L1 code audit at 3ddc30b: DONE (read-only).** See `audit_3ddc30b.json` for
   SHA-256 values and line ranges.

## Code audit findings

- `spectre_fo_enter` (hybrid.f90 lines 78-149) selects these once, from the
  entry state: layer width, h_entry=0.1 ro0/bmod, refine length, gmin, E0=|v|²/2,
  m_ref (radial_scale) and h_tau. The gyrophase and the `gauss2_encode`
  coordinates are also entry choices, and the encode depends on the selected
  width. All of this is the selection sigma(z0).
- No parameter changes inside a segment. None of these fields is assigned
  outside `enter`. There is no step halving and no retry after a Newton failure:
  the segment ends with STEP_LIMIT and recovery goes back to GC. The step h is
  never clipped. When a step overshoots the output clock, the excess is carried
  in `fo_time_ahead`. A segment survives across microsteps through
  `fo%active`. ro0 is global, the Newton tolerance is fixed at 1e-12 (and
  `fo_tolerance` is unused by Gauss2), and hfac=1 in production.
- Coverage: the in-segment steps are covered by the fixed-lambda leaf theorem,
  subject to the L4 caveats (Newton inexactness and branch regularity). These
  are not covered, because they are state-dependent switches: entry selection,
  the exit stopping rule (`check_exit`, which makes N state-dependent), the
  FO->GC map with recanonicalization, and re-entry.
- Minor: `gc_check_time` compares the physical clock against 0.25*h_tau, which
  is in tau units. This only changes how often exits are checked, not the map.

## Canonical entry choices (not proved)

Quantizing h, m_ref, E0 and the layer parameters to a fixed lattice makes sigma
locally constant away from the lattice switching set. Then
`locally_constant_selection` gives symplecticity there, but the lattice
boundaries remain switches. E0 only enters K = t(H-E0) as a shift. The
symplecticity of Phi_lam does not depend on its value, but selecting E0 = H(z0)
from the state is still selection.
