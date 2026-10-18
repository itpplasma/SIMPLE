import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Analysis.Calculus.FDeriv.Symmetric
import Mathlib.LinearAlgebra.BilinearMap
import Mathlib.Tactic

/- The tangent-lift quadratic invariant argument, with no preservation axiom.
   All stages and parameters here are fixed for the differentiation. -/
noncomputable section
namespace Gauss2
variable {V : Type*} [AddCommGroup V] [Module ℝ V]

def a11 : ℝ := 1/4
def a12 : ℝ := 1/4 - Real.sqrt 3 / 6
def a21 : ℝ := 1/4 + Real.sqrt 3 / 6
def a22 : ℝ := 1/4
def b : ℝ := 1/2

theorem butcher :
    b*a11 + b*a11 - b*b = 0 ∧
    b*a12 + b*a21 - b*b = 0 ∧
    b*a22 + b*a22 - b*b = 0 := by
  dsimp [b, a11, a12, a21, a22]
  constructor
  · ring
  constructor <;> ring

/- omega(u,v) is a quadratic invariant on the doubled tangent space V×V.
   This identity proves its preservation rather than assuming that RK
   schemes preserve it. The four Butcher defects appear explicitly. -/
theorem quadratic_defect (ω : V →ₗ[ℝ] V →ₗ[ℝ] ℝ)
    (u v k₁ k₂ l₁ l₂ : V) (h b₁ b₂ c₁₁ c₁₂ c₂₁ c₂₂ : ℝ) :
    ω (u + h • (b₁ • k₁ + b₂ • k₂))
      (v + h • (b₁ • l₁ + b₂ • l₂)) - ω u v =
    h*b₁*(ω k₁ (v+h•(c₁₁•l₁+c₁₂•l₂)) +
              ω (u+h•(c₁₁•k₁+c₁₂•k₂)) l₁) +
    h*b₂*(ω k₂ (v+h•(c₂₁•l₁+c₂₂•l₂)) +
              ω (u+h•(c₂₁•k₁+c₂₂•k₂)) l₂) -
    h^2*((b₁*c₁₁+b₁*c₁₁-b₁*b₁)*ω k₁ l₁ +
         (b₁*c₁₂+b₂*c₂₁-b₁*b₂)*ω k₁ l₂ +
         (b₂*c₂₁+b₁*c₁₂-b₂*b₁)*ω k₂ l₁ +
         (b₂*c₂₂+b₂*c₂₂-b₂*b₂)*ω k₂ l₂) := by
  simp only [map_add, map_smul, LinearMap.add_apply, LinearMap.smul_apply,
    smul_eq_mul]
  ring

theorem symplectic_tangent (ω : V →ₗ[ℝ] V →ₗ[ℝ] ℝ)
    (u v k₁ k₂ l₁ l₂ : V) (h : ℝ)
    (stage₁ : ω k₁ (v+h•(a11•l₁+a12•l₂)) +
              ω (u+h•(a11•k₁+a12•k₂)) l₁ = 0)
    (stage₂ : ω k₂ (v+h•(a21•l₁+a22•l₂)) +
              ω (u+h•(a21•k₁+a22•k₂)) l₂ = 0) :
    ω (u+h•(b•k₁+b•k₂)) (v+h•(b•l₁+b•l₂)) = ω u v := by
  have hd := quadratic_defect ω u v k₁ k₂ l₁ l₂ h b b a11 a12 a21 a22
  obtain ⟨h11,h12,h22⟩ := butcher
  have h21 : b*a21+b*a12-b*b=0 := by linarith
  rw [stage₁,stage₂,h11,h12,h21,h22] at hd
  linarith

/- C² gives the infinitesimal Hamiltonian pairing cancellation. omega is
   the constant canonical symplectic form, expressed as a continuous bilinear
   map; the Hamiltonian identity fixes the sign convention. -/
section Calculus
variable {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]

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
/- Exact stages on a differentiable branch. The explicit hY premises are
   where regularity of the implicit solve enters; no symplecticity premise. -/
theorem gauss2_symplectic
    (ω : E →L[ℝ] E →L[ℝ] ℝ) (H : E → ℝ) (F : E → E)
    (Y₁ Y₂ : E → E) (x : E) (h : ℝ)
    (S₁ S₂ D₁ D₂ : E →L[ℝ] E)
    (hY₁ : HasFDerivAt Y₁ S₁ x) (hY₂ : HasFDerivAt Y₂ S₂ x)
    (hF₁ : HasFDerivAt F D₁ (Y₁ x)) (hF₂ : HasFDerivAt F D₂ (Y₂ x))
    (hH₁ : ContDiffAt ℝ 2 H (Y₁ x)) (hH₂ : ContDiffAt ℝ 2 H (Y₂ x))
    (ham : ∀ y, ω (F y) = fderiv ℝ H y)
    (skew : ∀ u v, ω u v = -ω v u)
    (stage₁ : ∀ y, Y₁ y = y+h•(a11•F (Y₁ y)+a12•F (Y₂ y)))
    (stage₂ : ∀ y, Y₂ y = y+h•(a21•F (Y₁ y)+a22•F (Y₂ y))) :
    ∃ M : E →L[ℝ] E,
      HasFDerivAt (fun y => y+h•(b•F (Y₁ y)+b•F (Y₂ y))) M x ∧
      ∀ u v, ω (M u) (M v) = ω u v := by
  have hc₁ := hF₁.comp x hY₁
  have hc₂ := hF₂.comp x hY₂
  have hr₁ := (hasFDerivAt_id x).add
    (((hc₁.const_smul a11).add (hc₂.const_smul a12)).const_smul h)
  have hr₂ := (hasFDerivAt_id x).add
    (((hc₁.const_smul a21).add (hc₂.const_smul a22)).const_smul h)
  have he₁ : Y₁ = (fun y => y+h•(a11•F (Y₁ y)+a12•F (Y₂ y))) := funext stage₁
  have he₂ : Y₂ = (fun y => y+h•(a21•F (Y₁ y)+a22•F (Y₂ y))) := funext stage₂
  have hd₁ : S₁ = (ContinuousLinearMap.id ℝ E) +
      h • (a11 • (D₁.comp S₁) + a12 • (D₂.comp S₂)) := by
    exact (he₁ ▸ hY₁).unique hr₁
  have hd₂ : S₂ = (ContinuousLinearMap.id ℝ E) +
      h • (a21 • (D₁.comp S₁) + a22 • (D₂.comp S₂)) := by
    exact (he₂ ▸ hY₂).unique hr₂
  let M := (ContinuousLinearMap.id ℝ E) +
    h • (b • (D₁.comp S₁) + b • (D₂.comp S₂))
  refine ⟨M, (hasFDerivAt_id x).add
    (((hc₁.const_smul b).add (hc₂.const_smul b)).const_smul h), ?_⟩
  intro u v
  let B : E →ₗ[ℝ] E →ₗ[ℝ] ℝ := {
    toFun := fun u => (ω u).toLinearMap
    map_add' := by intro u v; ext w; simp
    map_smul' := by intro c u; ext w; simp }
  have hu₁ : S₁ u = u+h•(a11•D₁ (S₁ u)+a12•D₂ (S₂ u)) := by
    have hh := congrArg (fun L : E →L[ℝ] E => L u) hd₁
    simpa using hh
  have hv₁ : S₁ v = v+h•(a11•D₁ (S₁ v)+a12•D₂ (S₂ v)) := by
    have hh := congrArg (fun L : E →L[ℝ] E => L v) hd₁
    simpa using hh
  have hu₂ : S₂ u = u+h•(a21•D₁ (S₁ u)+a22•D₂ (S₂ u)) := by
    have hh := congrArg (fun L : E →L[ℝ] E => L u) hd₂
    simpa using hh
  have hv₂ : S₂ v = v+h•(a21•D₁ (S₁ v)+a22•D₂ (S₂ v)) := by
    have hh := congrArg (fun L : E →L[ℝ] E => L v) hd₂
    simpa using hh
  have h₁ := hamiltonian_linearization ω H F (Y₁ x) D₁ hH₁ hF₁ ham skew (S₁ u) (S₁ v)
  have h₂ := hamiltonian_linearization ω H F (Y₂ x) D₂ hH₂ hF₂ ham skew (S₂ u) (S₂ v)
  have ht₁ : B (D₁ (S₁ u)) (v+h•(a11•D₁ (S₁ v)+a12•D₂ (S₂ v))) +
      B (u+h•(a11•D₁ (S₁ u)+a12•D₂ (S₂ u))) (D₁ (S₁ v)) = 0 := by
    change ω (D₁ (S₁ u)) _ + ω _ (D₁ (S₁ v)) = 0
    rw [← hu₁, ← hv₁]
    exact h₁
  have ht₂ : B (D₂ (S₂ u)) (v+h•(a21•D₁ (S₁ v)+a22•D₂ (S₂ v))) +
      B (u+h•(a21•D₁ (S₁ u)+a22•D₂ (S₂ u))) (D₂ (S₂ v)) = 0 := by
    change ω (D₂ (S₂ u)) _ + ω _ (D₂ (S₂ v)) = 0
    rw [← hu₂, ← hv₂]
    exact h₂
  have hh := symplectic_tangent B u v (D₁ (S₁ u)) (D₂ (S₂ u))
    (D₁ (S₁ v)) (D₂ (S₂ v)) h ht₁ ht₂
  simpa [M, B] using hh

#print axioms gauss2_symplectic
end Calculus

#print axioms hamiltonian_linearization
#print axioms butcher
#print axioms quadratic_defect
#print axioms symplectic_tangent
end Gauss2
