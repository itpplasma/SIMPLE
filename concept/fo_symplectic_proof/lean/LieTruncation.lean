import Mathlib.Algebra.MvPolynomial.Basic
import Mathlib.Algebra.Polynomial.Eval.Defs
import Mathlib.Topology.Algebra.Polynomial
import Mathlib.Topology.Order.Compact
import Mathlib.Analysis.Normed.Group.Basic
import Mathlib.Analysis.Normed.Module.Basic
import Mathlib.Tactic

/- Truncation order of a pulled-back symplectic form (item D, GC <-> FO switch).

   What is proved (no sorry, no axioms beyond Lean/Mathlib's three):

   * `eval_cong`: composing a polynomial (any `MvPolynomial`) with two maps
     that agree modulo eps^k gives values that agree modulo eps^k.  This is
     the algebraic content of "pullback/composition of a truncated
     near-identity map agrees with the exact series through the truncation
     order" in a polynomial-in-eps model (the exact series is represented
     by any higher-order polynomial completion).
   * `pullback_cong`: the eps-scaled particle 2-form
        eps*Omega_p = F(x) dx^dx + eps * dv^dx
     pulled back by Jacobians whose position rows agree modulo eps^(N+2)
     and whose velocity rows agree modulo eps^(N+1) (and positions agree
     modulo eps^(N+2)) gives eps-scaled forms that agree modulo eps^(N+2).
     This is exactly the selector ordering of the GC prototype
     (position one order beyond velocity).
   * `defect_order_poly`: in R[eps], if eps*(pulled-back form - target)
     is divisible by eps^(N+2), the unscaled defect is eps^(N+1) times a
     polynomial and is bounded by C |eps|^(N+1) on every |eps| <= eps0.
   * `defect_bound`: the normed-space version for a non-polynomial family
     Omega(eps) = target(eps) + eps^(N+1) R(eps) with ||R|| <= C.

   What is NOT proved here: that the guiding-centre Lie transform makes the
   exact (untruncated) pullback equal to the GC form modulo an exact
   differential (that is the Lie-transform derivation, an input), that
   derivatives of truncation remainders inherit the eps order (a hypothesis
   on Jacobian rows here), or any bound for non-polynomial fields beyond the
   stated remainder hypothesis. -/

namespace LieTruncation

section Cong
variable {S : Type*} [CommRing S]

/-- `a` and `b` agree modulo `eps^k`. -/
def Cong (eps : S) (k : ℕ) (a b : S) : Prop := eps ^ k ∣ a - b

theorem cong_refl (eps : S) (k : ℕ) (a : S) : Cong eps k a a := by
  simp [Cong]

theorem cong_add {eps : S} {k : ℕ} {a b c d : S}
    (h1 : Cong eps k a b) (h2 : Cong eps k c d) : Cong eps k (a + c) (b + d) := by
  unfold Cong at *
  have : a + c - (b + d) = (a - b) + (c - d) := by ring
  rw [this]
  exact dvd_add h1 h2

theorem cong_sub {eps : S} {k : ℕ} {a b c d : S}
    (h1 : Cong eps k a b) (h2 : Cong eps k c d) : Cong eps k (a - c) (b - d) := by
  unfold Cong at *
  have : a - c - (b - d) = (a - b) - (c - d) := by ring
  rw [this]
  exact dvd_sub h1 h2

theorem cong_mul {eps : S} {k : ℕ} {a b c d : S}
    (h1 : Cong eps k a b) (h2 : Cong eps k c d) : Cong eps k (a * c) (b * d) := by
  unfold Cong at *
  have : a * c - b * d = (a - b) * c + b * (c - d) := by ring
  rw [this]
  exact dvd_add (dvd_mul_of_dvd_left h1 c) (dvd_mul_of_dvd_right h2 b)

theorem cong_mono {eps : S} {k l : ℕ} {a b : S} (hkl : k ≤ l)
    (h : Cong eps l a b) : Cong eps k a b :=
  dvd_trans (pow_dvd_pow eps hkl) h

/-- Multiplying by eps raises the order by one. -/
theorem cong_eps_mul {eps : S} {k : ℕ} {a b : S}
    (h : Cong eps k a b) : Cong eps (k + 1) (eps * a) (eps * b) := by
  unfold Cong at *
  obtain ⟨c, hc⟩ := h
  refine ⟨c, ?_⟩
  rw [← mul_sub, hc]
  ring

theorem cong_sum {ι : Type*} (s : Finset ι) {eps : S} {k : ℕ} (f g : ι → S)
    (h : ∀ i ∈ s, Cong eps k (f i) (g i)) :
    Cong eps k (∑ i ∈ s, f i) (∑ i ∈ s, g i) := by
  unfold Cong at *
  rw [← Finset.sum_sub_distrib]
  exact Finset.dvd_sum h

/-- Composition with a polynomial preserves agreement modulo eps^k. -/
theorem eval_cong {σ : Type*} {eps : S} {k : ℕ} (F : MvPolynomial σ S)
    (x y : σ → S) (h : ∀ i, Cong eps k (x i) (y i)) :
    Cong eps k (MvPolynomial.eval x F) (MvPolynomial.eval y F) := by
  induction F using MvPolynomial.induction_on with
  | h_C a => simp [Cong]
  | h_add p q hp hq =>
      simpa only [map_add] using cong_add hp hq
  | h_X p i hp =>
      simpa only [map_mul, MvPolynomial.eval_X] using cong_mul hp (h i)

/-- eps-scaled pullback of the particle form
    eps*Omega_p(a,b) = sum_kl a_k F_kl(x) b_l + eps*sum_k (av_k bx_k - ax_k bv_k)
    by the Jacobian with position rows `Jx` and velocity rows `Jv`. -/
def pullback {P Z : Type*} [Fintype P] (eps : S) (F : P → P → MvPolynomial P S)
    (x : P → S) (Jx Jv : P → Z → S) (i j : Z) : S :=
  (∑ k, ∑ l, Jx k i * MvPolynomial.eval x (F k l) * Jx l j) +
    eps * ∑ k, (Jv k i * Jx k j - Jx k i * Jv k j)

/-- Truncation lemma: position data agreeing modulo eps^(N+2) and velocity
    rows modulo eps^(N+1) give eps-scaled pulled-back forms agreeing modulo
    eps^(N+2), i.e. unscaled forms agreeing through order eps^N. -/
theorem pullback_cong {P Z : Type*} [Fintype P] (eps : S) (N : ℕ)
    (F : P → P → MvPolynomial P S) (x x' : P → S) (Jx Jx' Jv Jv' : P → Z → S)
    (hx : ∀ k, Cong eps (N + 2) (x k) (x' k))
    (hJx : ∀ k i, Cong eps (N + 2) (Jx k i) (Jx' k i))
    (hJv : ∀ k i, Cong eps (N + 1) (Jv k i) (Jv' k i)) (i j : Z) :
    Cong eps (N + 2) (pullback eps F x Jx Jv i j) (pullback eps F x' Jx' Jv' i j) := by
  unfold pullback
  apply cong_add
  · apply cong_sum
    intro k _
    apply cong_sum
    intro l _
    exact cong_mul (cong_mul (hJx k i) (eval_cong (F k l) x x' hx)) (hJx l j)
  · apply cong_eps_mul
    apply cong_sum
    intro k _
    have hxi := cong_mono (Nat.le_succ (N + 1)) (hJx k i)
    have hxj := cong_mono (Nat.le_succ (N + 1)) (hJx k j)
    exact cong_sub (cong_mul (hJv k i) hxj) (cong_mul hxi (hJv k j))

/-- Combining with an exact identity for the reference (untruncated) map:
    if the reference pullback equals the eps-scaled target modulo eps^(N+2),
    so does the truncated one. -/
theorem truncated_matches_target {P Z : Type*} [Fintype P] (eps : S) (N : ℕ)
    (F : P → P → MvPolynomial P S) (x x' : P → S) (Jx Jx' Jv Jv' : P → Z → S)
    (target : Z → Z → S)
    (hx : ∀ k, Cong eps (N + 2) (x k) (x' k))
    (hJx : ∀ k i, Cong eps (N + 2) (Jx k i) (Jx' k i))
    (hJv : ∀ k i, Cong eps (N + 1) (Jv k i) (Jv' k i))
    (hexact : ∀ i j, Cong eps (N + 2) (pullback eps F x Jx Jv i j) (target i j))
    (i j : Z) :
    Cong eps (N + 2) (pullback eps F x' Jx' Jv' i j) (target i j) := by
  have h := pullback_cong eps N F x x' Jx Jx' Jv Jv' hx hJx hJv i j
  unfold Cong at *
  have : pullback eps F x' Jx' Jv' i j - target i j =
      (pullback eps F x Jx Jv i j - target i j) -
        (pullback eps F x Jx Jv i j - pullback eps F x' Jx' Jv' i j) := by ring
  rw [this]
  exact dvd_sub (hexact i j) h

end Cong

section Bounds
open Polynomial

/-- In R[eps]: an eps-scaled defect divisible by eps^(N+2) is eps times an
    unscaled defect of the form eps^(N+1) c, bounded by C |eps|^(N+1) on any
    interval |eps| <= eps0. -/
theorem defect_order_poly (N : ℕ) (Δ : ℝ[X]) (h : X ^ (N + 2) ∣ Δ) (eps0 : ℝ) :
    ∃ δ : ℝ[X], Δ = X * δ ∧ ∃ C : ℝ, ∀ e : ℝ, |e| ≤ eps0 →
      |δ.eval e| ≤ C * |e| ^ (N + 1) := by
  obtain ⟨c, hc⟩ := h
  refine ⟨X ^ (N + 1) * c, ?_, ?_⟩
  · rw [hc]; ring
  · obtain ⟨C, hC⟩ := (isCompact_Icc (a := -eps0) (b := eps0)).exists_bound_of_continuousOn
      (c.continuous.continuousOn)
    refine ⟨C, fun e he => ?_⟩
    have hmem : e ∈ Set.Icc (-eps0) eps0 := abs_le.mp he
    have hb := hC e hmem
    rw [Real.norm_eq_abs] at hb
    rw [eval_mul, eval_pow, eval_X, abs_mul, abs_pow, mul_comm]
    exact mul_le_mul_of_nonneg_right hb (pow_nonneg (abs_nonneg e) _)

/-- Non-polynomial family: Omega(eps) = target(eps) + eps^(N+1) R(eps) with
    ||R(eps)|| <= C on |eps| <= eps0 gives a defect <= C |eps|^(N+1). -/
theorem defect_bound {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    (N : ℕ) (Ω target R : ℝ → E) (C eps0 : ℝ)
    (hsplit : ∀ e, Ω e = target e + e ^ (N + 1) • R e)
    (hR : ∀ e, |e| ≤ eps0 → ‖R e‖ ≤ C) :
    ∀ e, |e| ≤ eps0 → ‖Ω e - target e‖ ≤ C * |e| ^ (N + 1) := by
  intro e he
  rw [hsplit e, add_sub_cancel_left, norm_smul, Real.norm_eq_abs, abs_pow, mul_comm]
  exact mul_le_mul_of_nonneg_right (hR e he) (pow_nonneg (abs_nonneg e) _)

end Bounds
end LieTruncation

#print axioms LieTruncation.eval_cong
#print axioms LieTruncation.pullback_cong
#print axioms LieTruncation.truncated_matches_target
#print axioms LieTruncation.defect_order_poly
#print axioms LieTruncation.defect_bound
