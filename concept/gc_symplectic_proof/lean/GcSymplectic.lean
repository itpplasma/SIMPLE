import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Analysis.SpecialFunctions.Sqrt
import Mathlib.Analysis.Calculus.FDeriv.Symmetric
import Mathlib.LinearAlgebra.BilinearMap
import Mathlib.Tactic

/-!
# Symplecticity of SIMPLE's guiding-centre integrators (exact-stage level)

All shipped GC variants are partitioned Runge-Kutta (PRK) maps of the canonical
system with q = (theta, phi), p = (p_theta, p_phi) and Hamiltonian
K(q, p) = H(R(q, p), q, p_phi), R the implicit radius:

* `EXPL_IMPL_EULER`: s = 1, a = 0, â = 1, b = 1 (evaluation at (q0, p1)).
* `IMPL_EXPL_EULER`: s = 1, a = 1, â = 0, b = 1 (evaluation at (q1, p0)).
* `MIDPOINT`, `GAUSS1..4`: RK, i.e. PRK with â = a (Gauss-Legendre tableaux).
* `LOBATTO3`: a = Lobatto IIIA on q, â = Lobatto IIIB on p, b = (1,4,1)/6.

One generic theorem (`prk_symplectic`) proves that every PRK whose tableau
satisfies b_i â_ij + b_j a_ji - b_i b_j = 0 (with b̂ = b) has a symplectic
derivative on an exact differentiable stage branch. The conditions are then
checked exactly for each tableau the code ships. No symplecticity, invariant
preservation or generating-function axiom is assumed.
-/

noncomputable section
open Finset BigOperators Filter Topology

namespace GcSymplectic

/-! ## Scalar core: the O(h²) Butcher defect cancels -/

theorem prk_scalar {s : ℕ} (a ah : Fin s → Fin s → ℝ) (b : Fin s → ℝ)
    (hc : ∀ i j, b i * ah i j + b j * a j i - b i * b j = 0)
    (X Y : Fin s → Fin s → ℝ) (α β : Fin s → ℝ) (h : ℝ)
    (stage : ∀ i, α i + β i + h * ∑ j, (a i j * Y i j + ah i j * X i j)
        + h * ∑ j, (a i j * X j i + ah i j * Y j i) = 0) :
    h * ∑ i, b i * α i + h * ∑ i, b i * β i
      + h ^ 2 * ∑ i, ∑ j, b i * b j * (X i j + Y i j) = 0 := by
  have e1 : ∀ i, b i * α i + b i * β i =
      -(h * (∑ j, b i * (a i j * Y i j + ah i j * X i j)
        + ∑ j, b i * (a i j * X j i + ah i j * Y j i))) := by
    intro i
    rw [← Finset.mul_sum, ← Finset.mul_sum]
    linear_combination b i * stage i
  have swap : ∑ i, ∑ j, b i * (a i j * X j i + ah i j * Y j i) =
      ∑ i, ∑ j, b j * (a j i * X i j + ah j i * Y i j) := Finset.sum_comm
  have zero : ∑ i, ∑ j, (b i * b j * (X i j + Y i j)
      - b i * (a i j * Y i j + ah i j * X i j)
      - b j * (a j i * X i j + ah j i * Y i j)) = 0 := by
    apply Finset.sum_eq_zero; intro i _
    apply Finset.sum_eq_zero; intro j _
    linear_combination (-X i j) * hc i j + (-Y i j) * hc j i
  have hsum : ∑ i, b i * α i + ∑ i, b i * β i =
      -(h * (∑ i, ∑ j, b i * (a i j * Y i j + ah i j * X i j)
        + ∑ i, ∑ j, b i * (a i j * X j i + ah i j * Y j i))) := by
    rw [← Finset.sum_add_distrib, Finset.sum_congr rfl (fun i _ => e1 i),
      Finset.sum_neg_distrib, ← Finset.mul_sum, Finset.sum_add_distrib]
  rw [swap] at hsum
  simp only [Finset.sum_sub_distrib] at zero
  linear_combination h * hsum + h ^ 2 * zero

/-! ## Tangent-level PRK identity on a Lagrangian splitting -/

section Algebra
variable {V : Type*} [AddCommGroup V] [Module ℝ V]

theorem bilinear_expand {s : ℕ} (ω : V →ₗ[ℝ] V →ₗ[ℝ] ℝ) (Pq Pp : V →ₗ[ℝ] V)
    (hsplit : ∀ u, Pq u + Pp u = u)
    (lagq : ∀ u v, ω (Pq u) (Pq v) = 0) (lagp : ∀ u v, ω (Pp u) (Pp v) = 0)
    (b : Fin s → ℝ) (u v : V) (k l : Fin s → V) (h : ℝ) :
    ω (u + h • ∑ i, b i • k i) (v + h • ∑ i, b i • l i) =
      ω u v + h * ∑ i, b i * ω (k i) v + h * ∑ i, b i * ω u (l i)
      + h ^ 2 * ∑ i, ∑ j, b i * b j *
          (ω (Pq (k i)) (Pp (l j)) + ω (Pp (k i)) (Pq (l j))) := by
  have hkl : ∀ i j, ω (k i) (l j) =
      ω (Pq (k i)) (Pp (l j)) + ω (Pp (k i)) (Pq (l j)) := by
    intro i j
    conv_lhs => rw [← hsplit (k i), ← hsplit (l j)]
    simp only [map_add, LinearMap.add_apply, lagq, lagp]
    ring
  simp only [map_add, map_smul, map_sum, LinearMap.add_apply, LinearMap.smul_apply,
    LinearMap.coeFn_sum, Finset.sum_apply, smul_eq_mul, hkl]
  simp only [mul_add, Finset.sum_add_distrib, Finset.mul_sum]
  have sw : ∀ F : Fin s → Fin s → ℝ,
      ∑ x : Fin s, ∑ i : Fin s, h * (b x * (h * (b i * F i x))) =
      ∑ x : Fin s, ∑ i : Fin s, h ^ 2 * (b x * b i * F x i) := by
    intro F
    rw [Finset.sum_comm]
    apply Finset.sum_congr rfl; intro x _
    apply Finset.sum_congr rfl; intro i _
    ring
  rw [sw (fun i x => ω (Pq (k i)) (Pp (l x))), sw (fun i x => ω (Pp (k i)) (Pq (l x)))]
  ring

theorem stage_expand {s : ℕ} (ω : V →ₗ[ℝ] V →ₗ[ℝ] ℝ) (Pq Pp : V →ₗ[ℝ] V)
    (hsplit : ∀ u, Pq u + Pp u = u)
    (lagq : ∀ u v, ω (Pq u) (Pq v) = 0) (lagp : ∀ u v, ω (Pp u) (Pp v) = 0)
    (a ah : Fin s → Fin s → ℝ) (u v : V) (k l : Fin s → V) (h : ℝ) (i : Fin s) :
    ω (k i) (v + h • ∑ j, (a i j • Pq (l j) + ah i j • Pp (l j))) +
      ω (u + h • ∑ j, (a i j • Pq (k j) + ah i j • Pp (k j))) (l i) =
    ω (k i) v + ω u (l i)
      + h * ∑ j, (a i j * ω (Pp (k i)) (Pq (l j)) + ah i j * ω (Pq (k i)) (Pp (l j)))
      + h * ∑ j, (a i j * ω (Pq (k j)) (Pp (l i)) + ah i j * ω (Pp (k j)) (Pq (l i))) := by
  have h1 : ∀ x y, ω x (Pq y) = ω (Pp x) (Pq y) := by
    intro x y
    conv_lhs => rw [← hsplit x]
    simp [lagq]
  have h2 : ∀ x y, ω x (Pp y) = ω (Pq x) (Pp y) := by
    intro x y
    conv_lhs => rw [← hsplit x]
    simp [lagp]
  have h3 : ∀ x y, ω (Pq x) y = ω (Pq x) (Pp y) := by
    intro x y
    conv_lhs => rw [← hsplit y]
    simp [lagq]
  have h4 : ∀ x y, ω (Pp x) y = ω (Pp x) (Pq y) := by
    intro x y
    conv_lhs => rw [← hsplit y]
    simp [lagp]
  simp only [map_add, map_smul, map_sum, LinearMap.add_apply, LinearMap.smul_apply,
    LinearMap.coeFn_sum, Finset.sum_apply, smul_eq_mul]
  simp only [h1 (k i), h2 (k i), h3 _ (l i), h4 _ (l i)]
  ring

/-- PRK tangent identity: with b̂ = b and the condition
b_i â_ij + b_j a_ji - b_i b_j = 0, the linearized stage equations force
preservation of ω. q-increments use `a`, p-increments use `â`. -/
theorem prk_tangent {s : ℕ} (ω : V →ₗ[ℝ] V →ₗ[ℝ] ℝ) (Pq Pp : V →ₗ[ℝ] V)
    (hsplit : ∀ u, Pq u + Pp u = u)
    (lagq : ∀ u v, ω (Pq u) (Pq v) = 0) (lagp : ∀ u v, ω (Pp u) (Pp v) = 0)
    (a ah : Fin s → Fin s → ℝ) (b : Fin s → ℝ)
    (hc : ∀ i j, b i * ah i j + b j * a j i - b i * b j = 0)
    (u v : V) (k l : Fin s → V) (h : ℝ)
    (stage : ∀ i, ω (k i) (v + h • ∑ j, (a i j • Pq (l j) + ah i j • Pp (l j))) +
      ω (u + h • ∑ j, (a i j • Pq (k j) + ah i j • Pp (k j))) (l i) = 0) :
    ω (u + h • ∑ i, b i • k i) (v + h • ∑ i, b i • l i) = ω u v := by
  rw [bilinear_expand ω Pq Pp hsplit lagq lagp b u v k l h]
  have hs := prk_scalar a ah b hc (fun i j => ω (Pq (k i)) (Pp (l j)))
    (fun i j => ω (Pp (k i)) (Pq (l j))) (fun i => ω (k i) v) (fun i => ω u (l i)) h
    (by
      intro i
      have := stage i
      rw [stage_expand ω Pq Pp hsplit lagq lagp a ah u v k l h i] at this
      linear_combination this)
  linear_combination hs

end Algebra

/-! ## Calculus: exact stage branches give a symplectic derivative -/

section Calculus
variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]

/-- C² Hamiltonian ⇒ the linearized Hamiltonian vector field is ω-skew. -/
theorem hamiltonian_linearization
    (ω : E →L[ℝ] E →L[ℝ] ℝ) (H : E → ℝ) (F : E → E)
    (x : E) (D : E →L[ℝ] E)
    (hH : ContDiffAt ℝ 2 H x) (hF : HasFDerivAt F D x)
    (ham : ∀ y, ω (F y) = fderiv ℝ H y)
    (skew : ∀ u v, ω u v = -ω v u) (u v : E) :
    ω (D u) v + ω u (D v) = 0 := by
  have hd : DifferentiableAt ℝ (fderiv ℝ H) x :=
    (hH.fderiv_right (m := 1) (by norm_num)).differentiableAt (by norm_num)
  have lhs := ω.hasFDerivAt.comp x hF
  have rhs := hd.hasFDerivAt
  have heq : (fun y => ω (F y)) = (fun y => fderiv ℝ H y) := funext ham
  change HasFDerivAt (fun y => ω (F y)) (ω.comp D) x at lhs
  rw [heq] at lhs
  have eqder := lhs.unique rhs
  have huv : ω (D u) v = fderiv ℝ (fderiv ℝ H) x u v := by
    have hh := congrArg (fun L : E →L[ℝ] E →L[ℝ] ℝ => L u v) eqder
    simpa using hh
  have hvu : ω (D v) u = fderiv ℝ (fderiv ℝ H) x v u := by
    have hh := congrArg (fun L : E →L[ℝ] E →L[ℝ] ℝ => L v u) eqder
    simpa using hh
  have symm := hH.isSymmSndFDerivAt (by norm_num) u v
  rw [skew u (D v), huv, hvu, symm]
  ring

/-- Main theorem. `Y i` is a differentiable branch of exact stage solutions on a
neighbourhood of `x` (this is where regularity of the implicit solve enters);
`Pq`, `Pp` are the complementary Lagrangian projections onto the q- and
p-directions. The resulting one-step map has a derivative preserving ω. -/
theorem prk_symplectic {s : ℕ}
    (ω : E →L[ℝ] E →L[ℝ] ℝ) (H : E → ℝ) (F : E → E) (Pq Pp : E →L[ℝ] E)
    (hsplit : ∀ u, Pq u + Pp u = u)
    (lagq : ∀ u v, ω (Pq u) (Pq v) = 0) (lagp : ∀ u v, ω (Pp u) (Pp v) = 0)
    (a ah : Fin s → Fin s → ℝ) (b : Fin s → ℝ)
    (hc : ∀ i j, b i * ah i j + b j * a j i - b i * b j = 0)
    (Y : Fin s → E → E) (x : E) (h : ℝ) (S D : Fin s → E →L[ℝ] E)
    (hY : ∀ i, HasFDerivAt (Y i) (S i) x) (hF : ∀ i, HasFDerivAt F (D i) (Y i x))
    (hH : ∀ i, ContDiffAt ℝ 2 H (Y i x))
    (ham : ∀ y, ω (F y) = fderiv ℝ H y)
    (skew : ∀ u v, ω u v = -ω v u)
    (stage : ∀ i, Y i =ᶠ[𝓝 x]
      fun y => y + h • ∑ j, (a i j • Pq (F (Y j y)) + ah i j • Pp (F (Y j y)))) :
    ∃ M : E →L[ℝ] E,
      HasFDerivAt (fun y => y + h • ∑ i, b i • F (Y i y)) M x ∧
      ∀ u v, ω (M u) (M v) = ω u v := by
  have hc' : ∀ i, HasFDerivAt (fun y => F (Y i y)) ((D i).comp (S i)) x :=
    fun i => (hF i).comp x (hY i)
  have hrhs : ∀ i, HasFDerivAt
      (fun y => y + h • ∑ j, (a i j • Pq (F (Y j y)) + ah i j • Pp (F (Y j y))))
      (ContinuousLinearMap.id ℝ E + h • ∑ j, (a i j • Pq.comp ((D j).comp (S j))
        + ah i j • Pp.comp ((D j).comp (S j)))) x := by
    intro i
    apply (hasFDerivAt_id x).add
    apply HasFDerivAt.const_smul
    apply HasFDerivAt.sum
    intro j _
    exact ((Pq.hasFDerivAt.comp x (hc' j)).const_smul (a i j)).add
      ((Pp.hasFDerivAt.comp x (hc' j)).const_smul (ah i j))
  have hS : ∀ i, S i = ContinuousLinearMap.id ℝ E + h • ∑ j, (a i j •
      Pq.comp ((D j).comp (S j)) + ah i j • Pp.comp ((D j).comp (S j))) :=
    fun i => (hY i).unique ((hrhs i).congr_of_eventuallyEq (stage i))
  have hSu : ∀ i w, S i w = w + h • ∑ j, (a i j • Pq (D j (S j w)) +
      ah i j • Pp (D j (S j w))) := by
    intro i w
    have hh := congrArg (fun L : E →L[ℝ] E => L w) (hS i)
    simpa using hh
  let M : E →L[ℝ] E := ContinuousLinearMap.id ℝ E +
    h • ∑ i, b i • (D i).comp (S i)
  refine ⟨M, ?_, ?_⟩
  · apply (hasFDerivAt_id x).add
    apply HasFDerivAt.const_smul
    apply HasFDerivAt.sum
    intro i _
    exact (hc' i).const_smul (b i)
  intro u v
  let B : E →ₗ[ℝ] E →ₗ[ℝ] ℝ := {
    toFun := fun u => (ω u).toLinearMap
    map_add' := by intro u v; ext w; simp
    map_smul' := by intro c u; ext w; simp }
  have hM : ∀ w, M w = w + h • ∑ i, b i • D i (S i w) := by
    intro w; simp [M]
  have key := prk_tangent B Pq.toLinearMap Pp.toLinearMap hsplit lagq lagp a ah b hc
    u v (fun i => D i (S i u)) (fun i => D i (S i v)) h
    (by
      intro i
      have hl := hamiltonian_linearization ω H F (Y i x) (D i) (hH i) (hF i) ham skew
        (S i u) (S i v)
      change ω (D i (S i u)) (v + h • ∑ j, (a i j • Pq (D j (S j v)) +
          ah i j • Pp (D j (S j v)))) +
        ω (u + h • ∑ j, (a i j • Pq (D j (S j u)) + ah i j • Pp (D j (S j u))))
          (D i (S i v)) = 0
      rw [← hSu i u, ← hSu i v]
      exact hl)
  rw [hM u, hM v]
  exact key

/-- Classical Runge-Kutta corollary (â = a): b_i a_ij + b_j a_ji - b_i b_j = 0
suffices on the canonical splitting. -/
theorem rk_symplectic {s : ℕ}
    (ω : E →L[ℝ] E →L[ℝ] ℝ) (H : E → ℝ) (F : E → E) (Pq Pp : E →L[ℝ] E)
    (hsplit : ∀ u, Pq u + Pp u = u)
    (lagq : ∀ u v, ω (Pq u) (Pq v) = 0) (lagp : ∀ u v, ω (Pp u) (Pp v) = 0)
    (a : Fin s → Fin s → ℝ) (b : Fin s → ℝ)
    (hc : ∀ i j, b i * a i j + b j * a j i - b i * b j = 0)
    (Y : Fin s → E → E) (x : E) (h : ℝ) (S D : Fin s → E →L[ℝ] E)
    (hY : ∀ i, HasFDerivAt (Y i) (S i) x) (hF : ∀ i, HasFDerivAt F (D i) (Y i x))
    (hH : ∀ i, ContDiffAt ℝ 2 H (Y i x))
    (ham : ∀ y, ω (F y) = fderiv ℝ H y)
    (skew : ∀ u v, ω u v = -ω v u)
    (stage : ∀ i, Y i =ᶠ[𝓝 x] fun y => y + h • ∑ j, a i j • F (Y j y)) :
    ∃ M : E →L[ℝ] E,
      HasFDerivAt (fun y => y + h • ∑ i, b i • F (Y i y)) M x ∧
      ∀ u v, ω (M u) (M v) = ω u v := by
  refine prk_symplectic ω H F Pq Pp hsplit lagq lagp a a b hc Y x h S D hY hF hH ham
    skew ?_
  intro i
  refine (stage i).trans (Eventually.of_forall fun y => ?_)
  simp only
  congr 2
  apply Finset.sum_congr rfl
  intro j _
  rw [← smul_add, hsplit]

end Calculus

/-! ## Midpoint: identical to Gauss1 -/

section Midpoint
variable {W : Type*} [AddCommGroup W] [Module ℝ W]

/-- SIMPLE's `MIDPOINT` solves z1 = z0 + h F((z0+z1)/2); `GAUSS1` solves
Y = z0 + (h/2) F(Y), z1 = z0 + h F(Y). The solution sets coincide. -/
theorem midpoint_iff_gauss1 (F : W → W) (z0 z1 : W) (h : ℝ) :
    (∃ Y, Y = z0 + (h * (1 / 2)) • F Y ∧ z1 = z0 + (h * 1) • F Y) ↔
      z1 = z0 + h • F ((1 / 2 : ℝ) • (z0 + z1)) := by
  constructor
  · rintro ⟨Y, hY, h1⟩
    have hm : (1 / 2 : ℝ) • (z0 + z1) = Y := by
      calc (1 / 2 : ℝ) • (z0 + z1) = (1 / 2 : ℝ) • (z0 + (z0 + (h * 1) • F Y)) := by
            rw [h1]
        _ = z0 + (h * (1 / 2)) • F Y := by module
        _ = Y := hY.symm
    rw [hm]
    conv_lhs => rw [h1]
    module
  · intro h1
    refine ⟨(1 / 2 : ℝ) • (z0 + z1), ?_, ?_⟩
    · conv_lhs => rw [h1]
      module
    · conv_lhs => rw [h1]
      module

/-- The code's midpoint residuals fvec(5) and fvec(1) (both multiplied by
dpth(1) ≠ 0) force p_theta at the midpoint stage to be the average. -/
theorem midpoint_pth_average (dpm G p0 pm p1 h : ℝ) (hd : dpm ≠ 0)
    (f5 : dpm * (pm - p0) + (h / 2) * G = 0)
    (f1 : dpm * (p1 - p0) + h * G = 0) :
    pm = (p0 + p1) / 2 := by
  have : dpm * (2 * pm - p0 - p1) = 0 := by linear_combination 2 * f5 - f1
  rcases mul_eq_zero.mp this with h0 | h0
  · exact absurd h0 hd
  · linarith

end Midpoint

/-! ## Exact tableaux shipped by SIMPLE and their symplecticity conditions -/

def cond {s : ℕ} (a ah : Fin s → Fin s → ℝ) (b : Fin s → ℝ) : Prop :=
  ∀ i j, b i * ah i j + b j * a j i - b i * b j = 0

-- EXPL_IMPL_EULER: stage (q0, p1).
def euler1_a : Fin 1 → Fin 1 → ℝ := ![![0]]
def euler1_ah : Fin 1 → Fin 1 → ℝ := ![![1]]
def euler1_b : Fin 1 → ℝ := ![1]
-- IMPL_EXPL_EULER: stage (q1, p0).
def euler2_a : Fin 1 → Fin 1 → ℝ := ![![1]]
def euler2_ah : Fin 1 → Fin 1 → ℝ := ![![0]]
def euler2_b : Fin 1 → ℝ := ![1]

theorem euler1_cond : cond euler1_a euler1_ah euler1_b := by
  intro i j; fin_cases i; fin_cases j; simp [euler1_a, euler1_ah, euler1_b]
theorem euler2_cond : cond euler2_a euler2_ah euler2_b := by
  intro i j; fin_cases i; fin_cases j; simp [euler2_a, euler2_ah, euler2_b]

-- GAUSS1 (= MIDPOINT).
def gauss1_a : Fin 1 → Fin 1 → ℝ := ![![1 / 2]]
def gauss1_b : Fin 1 → ℝ := ![1]
-- GAUSS2.
def gauss2_a : Fin 2 → Fin 2 → ℝ :=
  ![![1 / 4, 1 / 4 - Real.sqrt 3 / 6],
    ![1 / 4 + Real.sqrt 3 / 6, 1 / 4]]
def gauss2_b : Fin 2 → ℝ := ![1 / 2, 1 / 2]
-- GAUSS3.
def gauss3_a : Fin 3 → Fin 3 → ℝ :=
  ![![5 / 36, 2 / 9 - Real.sqrt 15 / 15, 5 / 36 - Real.sqrt 15 / 30],
    ![5 / 36 + Real.sqrt 15 / 24, 2 / 9, 5 / 36 - Real.sqrt 15 / 24],
    ![5 / 36 + Real.sqrt 15 / 30, 2 / 9 + Real.sqrt 15 / 15, 5 / 36]]
def gauss3_b : Fin 3 → ℝ := ![5 / 18, 4 / 9, 5 / 18]
-- GAUSS4 (Butcher's closed form).
def g4t : ℝ := Real.sqrt 30
def g4u : ℝ := Real.sqrt ((15 + 2 * Real.sqrt 30) / 35)
def g4v : ℝ := Real.sqrt ((15 - 2 * Real.sqrt 30) / 35)
def g4w1 : ℝ := 1 / 8 - g4t / 144
def g4w1p : ℝ := 1 / 8 + g4t / 144
def g4w2 : ℝ := g4u / 2
def g4w2p : ℝ := g4v / 2
def g4w3 : ℝ := g4w2 * (1 / 6 + g4t / 24)
def g4w3p : ℝ := g4w2p * (1 / 6 - g4t / 24)
def g4w4 : ℝ := g4w2 * (1 / 21 + 5 * g4t / 168)
def g4w4p : ℝ := g4w2p * (1 / 21 - 5 * g4t / 168)
def g4w5 : ℝ := g4w2 - 2 * g4w3
def g4w5p : ℝ := g4w2p - 2 * g4w3p
def gauss4_a : Fin 4 → Fin 4 → ℝ :=
  ![![g4w1, g4w1p - g4w3 + g4w4p, g4w1p - g4w3 - g4w4p, g4w1 - g4w5],
    ![g4w1 - g4w3p + g4w4, g4w1p, g4w1p - g4w5p, g4w1 - g4w3p - g4w4],
    ![g4w1 + g4w3p + g4w4, g4w1p + g4w5p, g4w1p, g4w1 + g4w3p - g4w4],
    ![g4w1 + g4w5, g4w1p + g4w3 + g4w4p, g4w1p + g4w3 - g4w4p, g4w1]]
def gauss4_b : Fin 4 → ℝ := ![2 * g4w1, 2 * g4w1p, 2 * g4w1p, 2 * g4w1]
-- LOBATTO3: IIIA on q, IIIB on p.
def lobatto3_a : Fin 3 → Fin 3 → ℝ :=
  ![![0, 0, 0], ![5 / 24, 1 / 3, -1 / 24], ![1 / 6, 2 / 3, 1 / 6]]
def lobatto3_ah : Fin 3 → Fin 3 → ℝ :=
  ![![1 / 6, -1 / 6, 0], ![1 / 6, 1 / 3, 0], ![1 / 6, 5 / 6, 0]]
def lobatto3_b : Fin 3 → ℝ := ![1 / 6, 2 / 3, 1 / 6]

theorem gauss1_cond : cond gauss1_a gauss1_a gauss1_b := by
  intro i j; fin_cases i; fin_cases j; simp [gauss1_a, gauss1_b]; norm_num

theorem gauss2_cond : cond gauss2_a gauss2_a gauss2_b := by
  intro i j
  fin_cases i <;> fin_cases j <;> simp [gauss2_a, gauss2_b] <;> ring

theorem gauss3_cond : cond gauss3_a gauss3_a gauss3_b := by
  intro i j
  fin_cases i <;> fin_cases j <;> simp [gauss3_a, gauss3_b] <;> ring

theorem gauss4_cond : cond gauss4_a gauss4_a gauss4_b := by
  have ht : g4t ^ 2 = 30 := Real.sq_sqrt (by norm_num)
  intro i j
  fin_cases i <;> fin_cases j <;>
    simp [gauss4_a, gauss4_b, g4w1, g4w1p, g4w2, g4w2p, g4w3, g4w3p, g4w4, g4w4p,
      g4w5, g4w5p, Matrix.vecHead, Matrix.vecTail]
  -- rows i = 0..3, columns j = 0..3; coefficients from exact division by t²-30
  · linear_combination (0 : ℝ) * ht
  · linear_combination (g4u / 2016 + g4v / 2016) * ht
  · linear_combination (g4u / 2016 - g4v / 2016) * ht
  · linear_combination (0 : ℝ) * ht
  · linear_combination (g4u / 2016 + g4v / 2016) * ht
  · linear_combination (0 : ℝ) * ht
  · linear_combination (0 : ℝ) * ht
  · linear_combination (-g4u / 2016 + g4v / 2016) * ht
  · linear_combination (g4u / 2016 - g4v / 2016) * ht
  · linear_combination (0 : ℝ) * ht
  · linear_combination (0 : ℝ) * ht
  · linear_combination (-g4u / 2016 - g4v / 2016) * ht
  · linear_combination (0 : ℝ) * ht
  · linear_combination (-g4u / 2016 + g4v / 2016) * ht
  · linear_combination (-g4u / 2016 - g4v / 2016) * ht
  · linear_combination (0 : ℝ) * ht

theorem lobatto3_cond : cond lobatto3_a lobatto3_ah lobatto3_b := by
  intro i j
  fin_cases i <;> fin_cases j <;> simp [lobatto3_a, lobatto3_ah, lobatto3_b] <;> norm_num

/-- Consistency: the Gauss tableaux are not degenerate (row sums are the nodes
and Σ b = 1), so the conditions above are about genuine order-2s methods. -/
theorem gauss_weights_sum :
    ∑ i, gauss1_b i = 1 ∧ ∑ i, gauss2_b i = 1 ∧ ∑ i, gauss3_b i = 1 ∧
      ∑ i, lobatto3_b i = 1 ∧ ∑ i, gauss4_b i = 1 := by
  refine ⟨?_, ?_, ?_, ?_, ?_⟩ <;>
    simp [Fin.sum_univ_succ, gauss1_b, gauss2_b, gauss3_b, lobatto3_b, gauss4_b,
      g4w1, g4w1p] <;> ring


/-! ## Instantiation: every shipped GC variant -/

/-- The exact-stage statement of `prk_symplectic` for a fixed tableau. -/
def SymplecticOnBranches {s : ℕ} (a ah : Fin s → Fin s → ℝ) (b : Fin s → ℝ) : Prop :=
  ∀ (E : Type) [NormedAddCommGroup E] [NormedSpace ℝ E]
    (ω : E →L[ℝ] E →L[ℝ] ℝ) (H : E → ℝ) (F : E → E) (Pq Pp : E →L[ℝ] E),
    (∀ u, Pq u + Pp u = u) → (∀ u v, ω (Pq u) (Pq v) = 0) →
    (∀ u v, ω (Pp u) (Pp v) = 0) →
    ∀ (Y : Fin s → E → E) (x : E) (h : ℝ) (S D : Fin s → E →L[ℝ] E),
    (∀ i, HasFDerivAt (Y i) (S i) x) → (∀ i, HasFDerivAt F (D i) (Y i x)) →
    (∀ i, ContDiffAt ℝ 2 H (Y i x)) → (∀ y, ω (F y) = fderiv ℝ H y) →
    (∀ u v, ω u v = -ω v u) →
    (∀ i, Y i =ᶠ[𝓝 x]
      fun y => y + h • ∑ j, (a i j • Pq (F (Y j y)) + ah i j • Pp (F (Y j y)))) →
    ∃ M : E →L[ℝ] E,
      HasFDerivAt (fun y => y + h • ∑ i, b i • F (Y i y)) M x ∧
      ∀ u v, ω (M u) (M v) = ω u v

theorem symplectic_of_cond {s : ℕ} {a ah : Fin s → Fin s → ℝ} {b : Fin s → ℝ}
    (hc : cond a ah b) : SymplecticOnBranches a ah b := by
  intro E _ _ ω H F Pq Pp hsplit lagq lagp Y x h S D hY hF hH ham skew stage
  exact prk_symplectic ω H F Pq Pp hsplit lagq lagp a ah b hc Y x h S D hY hF hH ham
    skew stage

theorem expl_impl_euler_symplectic : SymplecticOnBranches euler1_a euler1_ah euler1_b :=
  symplectic_of_cond euler1_cond
theorem impl_expl_euler_symplectic : SymplecticOnBranches euler2_a euler2_ah euler2_b :=
  symplectic_of_cond euler2_cond
theorem gauss1_symplectic : SymplecticOnBranches gauss1_a gauss1_a gauss1_b :=
  symplectic_of_cond gauss1_cond
theorem gauss2_symplectic : SymplecticOnBranches gauss2_a gauss2_a gauss2_b :=
  symplectic_of_cond gauss2_cond
theorem gauss3_symplectic : SymplecticOnBranches gauss3_a gauss3_a gauss3_b :=
  symplectic_of_cond gauss3_cond
theorem gauss4_symplectic : SymplecticOnBranches gauss4_a gauss4_a gauss4_b :=
  symplectic_of_cond gauss4_cond
theorem lobatto3_symplectic : SymplecticOnBranches lobatto3_a lobatto3_ah lobatto3_b :=
  symplectic_of_cond lobatto3_cond

#print axioms prk_scalar
#print axioms prk_tangent
#print axioms hamiltonian_linearization
#print axioms prk_symplectic
#print axioms rk_symplectic
#print axioms midpoint_iff_gauss1
#print axioms midpoint_pth_average
#print axioms euler1_cond
#print axioms euler2_cond
#print axioms gauss1_cond
#print axioms gauss2_cond
#print axioms gauss3_cond
#print axioms gauss4_cond
#print axioms lobatto3_cond
#print axioms gauss_weights_sum
#print axioms expl_impl_euler_symplectic
#print axioms impl_expl_euler_symplectic
#print axioms gauss1_symplectic
#print axioms gauss2_symplectic
#print axioms gauss3_symplectic
#print axioms gauss4_symplectic
#print axioms lobatto3_symplectic

end GcSymplectic
