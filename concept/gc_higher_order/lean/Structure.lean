/- Abstract cancellation of invariant pairings under a flow-reset-flow sandwich.
   No axioms beyond the stated pairing-preservation and surjectivity hypotheses.
   To apply to differential geometry, V is a tangent bundle with basepoints
   encoded, and these maps are tangent lifts. This file does not establish
   differentiability, a GC expansion, energy conservation, or error bounds. -/

def Preserves {V R : Type} (ω : V → V → R) (f : V → V) : Prop :=
  ∀ v w, ω (f v) (f w) = ω v w

theorem preserves_comp {V R : Type} (ω : V → V → R)
    (f g : V → V) (hf : Preserves ω f) (hg : Preserves ω g) :
    Preserves ω (fun v => f (g v)) := by
  intro v w
  exact Eq.trans (hf (g v) (g w)) (hg v w)

theorem sandwich_iff {V R : Type} (ω : V → V → R)
    (after reset before : V → V)
    (ha : Preserves ω after) (hb : Preserves ω before)
    (hs : ∀ v, ∃ w, before w = v) :
    Preserves ω (fun v => after (reset (before v))) ↔ Preserves ω reset := by
  constructor
  · intro h v w
    obtain ⟨x, hx⟩ := hs v
    obtain ⟨y, hy⟩ := hs w
    rw [← hx, ← hy]
    exact Eq.trans (Eq.symm (ha (reset (before x)) (reset (before y))))
      (Eq.trans (h x y) (Eq.symm (hb x y)))
  · intro hr
    exact preserves_comp ω after (fun v => reset (before v)) ha
      (preserves_comp ω reset before hr hb)

#print axioms sandwich_iff
