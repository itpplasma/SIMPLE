import Gauss2

/- Entry-frozen parameter selection for the exact-stage Gauss2 map.

   1. Fixed-parameter composition: N exact-stage Gauss2 steps of one fixed
      Hamiltonian K_lam compose to a symplectic map (`frozen_composition`).
   2. Leaf/Poisson statement: (z,lam) ↦ (Phi_lam z, lam) keeps every leaf and
      is a Poisson map for the Casimir-extended bivector iff the leaf
      derivative preserves the leaf bivector; D_lam Phi never enters
      (`leaf_poisson_iff`, `leaf_poisson_of_symplectic`).
   3. Selection obstruction: Psi z = Phi (sigma z, z) has derivative
      D_lam Phi ∘ D sigma + D_z Phi (`selection_hasFDerivAt`) and the exact
      symplectic defect consists of the cross terms (`selection_defect`).
      Locally constant sigma gives a symplectic Psi
      (`locally_constant_selection`); a step switch is discontinuous
      (`switch_not_continuous`, `switch_not_symplectic`).
   4. Counterexample: H = p^2/2, actual Gauss2 coefficients, state-selected
      h(q): the symplectic pairing is 1 + h'(q) p (`selected_pairing`), and
      h = h0 + q, p = 1 is not symplectic (`selected_not_symplectic`).
   No sorry, no custom axiom. -/
noncomputable section
namespace EntryFrozen
open Gauss2

section General
variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
variable {P : Type*} [NormedAddCommGroup P] [NormedSpace ℝ P]

/-- Phi is symplectic at x: it has a derivative there that preserves omega. -/
def SymplAt (ω : E →L[ℝ] E →L[ℝ] ℝ) (Φ : E → E) (x : E) : Prop :=
  ∃ M : E →L[ℝ] E, HasFDerivAt Φ M x ∧ ∀ u v, ω (M u) (M v) = ω u v

theorem symplAt_id (ω : E →L[ℝ] E →L[ℝ] ℝ) (x : E) : SymplAt ω id x :=
  ⟨ContinuousLinearMap.id ℝ E, hasFDerivAt_id x, fun _ _ => rfl⟩

theorem symplAt_comp {ω : E →L[ℝ] E →L[ℝ] ℝ} {Φ Ψ : E → E} {x : E}
    (hΦ : SymplAt ω Φ x) (hΨ : SymplAt ω Ψ (Φ x)) : SymplAt ω (Ψ ∘ Φ) x := by
  obtain ⟨M, hM, hMω⟩ := hΦ
  obtain ⟨N, hN, hNω⟩ := hΨ
  refine ⟨N.comp M, hN.comp x hM, fun u v => ?_⟩
  rw [ContinuousLinearMap.comp_apply, ContinuousLinearMap.comp_apply, hNω, hMω]

/-- compose Φ N = Φ (N-1) ∘ ... ∘ Φ 0. -/
def compose (Φ : ℕ → E → E) : ℕ → E → E
  | 0 => id
  | n + 1 => Φ n ∘ compose Φ n

theorem compose_symplectic (ω : E →L[ℝ] E →L[ℝ] ℝ) (Φ : ℕ → E → E) (x : E) :
    ∀ N : ℕ, (∀ k, k < N → SymplAt ω (Φ k) (compose Φ k x)) →
      SymplAt ω (compose Φ N) x := by
  intro N
  induction N with
  | zero => intro _; exact symplAt_id ω x
  | succ n ih =>
    intro hstep
    have h1 := ih (fun k hk => hstep k (by omega))
    have h2 := hstep n (by omega)
    exact symplAt_comp h1 h2

/-- One exact-stage Gauss2 step of the vector field F with stages Y₁, Y₂. -/
def gaussStep (F : E → E) (Y₁ Y₂ : E → E) (h : ℝ) : E → E :=
  fun y => y + h • (b • F (Y₁ y) + b • F (Y₂ y))

/-- Claim 1. One fixed parameter lam fixes K = Hf lam, its vector field
    Ff lam and the step hf lam for all N steps. Stage branches may differ per
    step (each must be an exact differentiable branch along the orbit). -/
theorem frozen_composition (ω : E →L[ℝ] E →L[ℝ] ℝ)
    (Hf : P → E → ℝ) (Ff : P → E → E) (hf : P → ℝ) (lam : P)
    (Y₁ Y₂ : ℕ → E → E) (x : E) (N : ℕ)
    (S₁ S₂ D₁ D₂ : ℕ → E →L[ℝ] E)
    (ham : ∀ y, ω (Ff lam y) = fderiv ℝ (Hf lam) y)
    (skew : ∀ u v, ω u v = -ω v u)
    (stage₁ : ∀ k y, Y₁ k y = y + hf lam • (a11 • Ff lam (Y₁ k y) + a12 • Ff lam (Y₂ k y)))
    (stage₂ : ∀ k y, Y₂ k y = y + hf lam • (a21 • Ff lam (Y₁ k y) + a22 • Ff lam (Y₂ k y)))
    (hreg : ∀ k, k < N →
      let xk := compose (fun j => gaussStep (Ff lam) (Y₁ j) (Y₂ j) (hf lam)) k x
      HasFDerivAt (Y₁ k) (S₁ k) xk ∧ HasFDerivAt (Y₂ k) (S₂ k) xk ∧
      HasFDerivAt (Ff lam) (D₁ k) (Y₁ k xk) ∧ HasFDerivAt (Ff lam) (D₂ k) (Y₂ k xk) ∧
      ContDiffAt ℝ 2 (Hf lam) (Y₁ k xk) ∧ ContDiffAt ℝ 2 (Hf lam) (Y₂ k xk)) :
    SymplAt ω (compose (fun j => gaussStep (Ff lam) (Y₁ j) (Y₂ j) (hf lam)) N) x := by
  apply compose_symplectic
  intro k hk
  obtain ⟨a1, a2, a3, a4, a5, a6⟩ := hreg k hk
  exact gauss2_symplectic ω (Hf lam) (Ff lam) (Y₁ k) (Y₂ k) _ (hf lam)
    (S₁ k) (S₂ k) (D₁ k) (D₂ k) a1 a2 a3 a4 a5 a6 ham skew (stage₁ k) (stage₂ k)

/-! Claim 2: leaves and the Casimir-extended Poisson structure. -/

/-- The extended map G(z, lam) = (Phi lam z, lam). -/
def extend (Φ : P → E → E) : E × P → E × P := fun w => (Φ w.2 w.1, w.2)

theorem extend_leaf (Φ : P → E → E) (w : E × P) : (extend Φ w).2 = w.2 := rfl

theorem extend_on_leaf (Φ : P → E → E) (lam : P) :
    (fun z => (extend Φ (z, lam)).1) = Φ lam := rfl

/-- Derivative of G from a joint derivative L of (z,lam) ↦ Phi lam z. -/
theorem extend_hasFDerivAt (Φ : P → E → E) (w : E × P) (L : E × P →L[ℝ] E)
    (hΦ : HasFDerivAt (fun w : E × P => Φ w.2 w.1) L w) :
    HasFDerivAt (extend Φ) (L.prod (ContinuousLinearMap.snd ℝ E P)) w :=
  hΦ.prod (hasFDerivAt_snd)

/-- Casimir-extended bivector: pairs only the E-parts of covectors. -/
def extPi (πE : (E →L[ℝ] ℝ) → (E →L[ℝ] ℝ) → ℝ) (α β : E × P →L[ℝ] ℝ) : ℝ :=
  πE (α.comp (ContinuousLinearMap.inl ℝ E P)) (β.comp (ContinuousLinearMap.inl ℝ E P))

theorem pullback_inl (L : E × P →L[ℝ] E) (α : E × P →L[ℝ] ℝ) :
    (α.comp (L.prod (ContinuousLinearMap.snd ℝ E P))).comp
        (ContinuousLinearMap.inl ℝ E P) =
      (α.comp (ContinuousLinearMap.inl ℝ E P)).comp
        (L.comp (ContinuousLinearMap.inl ℝ E P)) := by
  ext u
  simp

/-- Claim 2 (Poisson form). G is a Poisson map of the Casimir-extended
    bivector iff A = D_z Phi preserves the leaf bivector. The parameter
    derivative L ∘ inr (= D_lam Phi) does not occur. -/
theorem leaf_poisson_iff (πE : (E →L[ℝ] ℝ) → (E →L[ℝ] ℝ) → ℝ)
    (L : E × P →L[ℝ] E) :
    (∀ α β : E × P →L[ℝ] ℝ,
        extPi πE (α.comp (L.prod (ContinuousLinearMap.snd ℝ E P)))
          (β.comp (L.prod (ContinuousLinearMap.snd ℝ E P))) = extPi πE α β) ↔
      (∀ α β : E →L[ℝ] ℝ,
        πE (α.comp (L.comp (ContinuousLinearMap.inl ℝ E P)))
          (β.comp (L.comp (ContinuousLinearMap.inl ℝ E P))) = πE α β) := by
  constructor
  · intro hG α β
    have hh := hG (α.comp (ContinuousLinearMap.fst ℝ E P))
      (β.comp (ContinuousLinearMap.fst ℝ E P))
    have e1 : ∀ γ : E →L[ℝ] ℝ,
        (γ.comp (ContinuousLinearMap.fst ℝ E P)).comp
          (ContinuousLinearMap.inl ℝ E P) = γ := by
      intro γ; ext u; simp
    unfold extPi at hh
    rw [pullback_inl, pullback_inl, e1, e1] at hh
    exact hh
  · intro hA α β
    unfold extPi
    rw [pullback_inl, pullback_inl]
    exact hA _ _

/-- The leaf bivector of a nondegenerate omega via a sharp map, and its
    preservation by a surjective omega-preserving A. -/
theorem leaf_poisson_of_symplectic (ω : E →L[ℝ] E →L[ℝ] ℝ)
    (sharp : (E →L[ℝ] ℝ) → E) (hsharp : ∀ β v, ω (sharp β) v = β v)
    (nondeg : ∀ x, (∀ v, ω x v = 0) → x = 0)
    (A : E →L[ℝ] E) (hA : ∀ u v, ω (A u) (A v) = ω u v)
    (hsurj : Function.Surjective A) (α β : E →L[ℝ] ℝ) :
    (α.comp A) (sharp (β.comp A)) = α (sharp β) := by
  have key : A (sharp (β.comp A)) = sharp β := by
    have hz := nondeg (A (sharp (β.comp A)) - sharp β) (by
      intro v
      obtain ⟨w, rfl⟩ := hsurj v
      rw [map_sub, ContinuousLinearMap.sub_apply, hA, hsharp, hsharp]
      simp)
    exact sub_eq_zero.mp hz
  rw [ContinuousLinearMap.comp_apply, key]

/-! Claim 3: state-dependent selection. -/

/-- Exact defect identity for M = C + A, C = D_lam Phi ∘ D sigma, A = D_z Phi. -/
theorem selection_defect (ω : E →L[ℝ] E →L[ℝ] ℝ) (A C : E →L[ℝ] E) (u v : E) :
    ω ((C + A) u) ((C + A) v) - ω u v =
      (ω (A u) (A v) - ω u v) + ω (A u) (C v) + ω (C u) (A v) + ω (C u) (C v) := by
  simp only [ContinuousLinearMap.add_apply, map_add]
  ring

/-- With symplectic A the defect is exactly the three cross terms. -/
theorem selection_defect_symplectic (ω : E →L[ℝ] E →L[ℝ] ℝ) (A C : E →L[ℝ] E)
    (hA : ∀ u v, ω (A u) (A v) = ω u v) (u v : E) :
    ω ((C + A) u) ((C + A) v) - ω u v =
      ω (A u) (C v) + ω (C u) (A v) + ω (C u) (C v) := by
  rw [selection_defect, hA]
  ring

/-- Locally constant selection (e.g. a lattice-quantized parameter away from
    its switching set): Psi agrees with a fixed-parameter map near x. -/
theorem locally_constant_selection (ω : E →L[ℝ] E →L[ℝ] ℝ) (Φ : P × E → E)
    (σ : E → P) (x : E) (c : P)
    (hloc : ∀ᶠ z in nhds x, σ z = c)
    (hfix : SymplAt ω (fun z => Φ (c, z)) x) :
    SymplAt ω (fun z => Φ (σ z, z)) x := by
  obtain ⟨M, hM, hMω⟩ := hfix
  refine ⟨M, hM.congr_of_eventuallyEq ?_, hMω⟩
  filter_upwards [hloc] with z hz
  rw [hz]

end General

#print axioms symplAt_comp
#print axioms compose_symplectic
#print axioms frozen_composition
#print axioms extend_leaf
#print axioms extend_on_leaf
#print axioms extend_hasFDerivAt
#print axioms pullback_inl
#print axioms leaf_poisson_iff
#print axioms leaf_poisson_of_symplectic
#print axioms selection_defect
#print axioms selection_defect_symplectic
#print axioms locally_constant_selection
end EntryFrozen
