/-
HalfGap.lean — Lemma 3 of work/DT3-THEOREM.md (half-gap identity V_I(h) = ½ sup (h − GCM h)) in a form that needs no
convex-minorant machinery (Step 8, DT6). Both directions of the proof are checked:

1. `fit_of_minorant` : if c is a convex minorant of h on s with gap h − c ≤ G, then c + G/2 is convex and within G/2 of h.
2. `minorant_of_fit` : if g is convex and within d of h on s, then g − d is a convex minorant of h with gap ≤ 2d.
Together: the best convex-fit error equals half the best minorant gap (the GCM is the best minorant), which is Lemma 3.
-/
import Mathlib

/-- **(≤) A convex minorant with gap G gives a convex fit within G/2.** -/
theorem fit_of_minorant {s : Set ℝ} {h c : ℝ → ℝ} {G : ℝ} (hc : ConvexOn ℝ s c)
    (hmin : ∀ x ∈ s, c x ≤ h x) (hgap : ∀ x ∈ s, h x - c x ≤ G) :
    ConvexOn ℝ s (fun x => c x + G / 2) ∧ ∀ x ∈ s, |h x - (c x + G / 2)| ≤ G / 2 := by
  refine ⟨hc.add_const (G / 2), fun x hx => ?_⟩
  rw [abs_le]
  constructor <;> linarith [hmin x hx, hgap x hx]

/-- **(≥) A convex fit within d gives a convex minorant with gap at most 2d.** -/
theorem minorant_of_fit {s : Set ℝ} {h g : ℝ → ℝ} {d : ℝ} (hg : ConvexOn ℝ s g)
    (hfit : ∀ x ∈ s, |h x - g x| ≤ d) :
    ConvexOn ℝ s (fun x => g x - d) ∧ (∀ x ∈ s, g x - d ≤ h x) ∧ ∀ x ∈ s, h x - (g x - d) ≤ 2 * d := by
  refine ⟨⟨hg.1, ?_⟩, fun x hx => ?_, fun x hx => ?_⟩
  · intro x hx y hy a b ha hb hab
    have hc := hg.2 hx hy ha hb hab
    simp only [smul_eq_mul] at *
    have hd : a * d + b * d = d := by rw [← add_mul, hab, one_mul]
    linarith
  · have := (abs_le.1 (hfit x hx)).1
    linarith
  · have := (abs_le.1 (hfit x hx)).2
    linarith

/-! NONVAC: h = x², c = 0 on [0, 1] (convex minorant, gap ≤ 1); the fit 1/2 is within 1/2. -/
example : ConvexOn ℝ (Set.Icc (0 : ℝ) 1) (fun x => (fun _ : ℝ => (0 : ℝ)) x + 1 / 2) ∧
    ∀ x ∈ Set.Icc (0 : ℝ) 1, |x ^ 2 - ((fun _ : ℝ => (0 : ℝ)) x + 1 / 2)| ≤ 1 / 2 :=
  fit_of_minorant (convexOn_const 0 (convex_Icc 0 1)) (fun x _ => sq_nonneg x)
    (fun x hx => by
      simp only [sub_zero]
      nlinarith [hx.1, hx.2])

#print axioms fit_of_minorant
#print axioms minorant_of_fit
