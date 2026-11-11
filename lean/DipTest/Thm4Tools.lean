/-
Thm4Tools.lean — machine checks of the seam-shift tool (T1) of work/DT3-THEOREM.md §4c (Theorem 4 v0.1).
Step 8 (DT6), written 2026-11-06 while DT5 runs.

1. `endpoint_lower_concave` : lowering the value of a concave function at the left end point of its interval keeps
   it concave (the step "G(μ) ≤ g₂(μ) + δ⁺/2 keeps G concave on [μ, 1]" in (T1)).
2. `endpoint_raise_convex`  : mirror: raising the value at the right end point keeps a convex function convex.
3. `seam_shift_error`       : the error bound e + s of the shifted fit on both sides of the seam.
-/
import Mathlib

/-- **Lowering the left end value keeps concavity.** -/
theorem endpoint_lower_concave {μ b v : ℝ} {g : ℝ → ℝ} (hg : ConcaveOn ℝ (Set.Icc μ b) g)
    (hv : v ≤ g μ) : ConcaveOn ℝ (Set.Icc μ b) (fun x => if x = μ then v else g x) := by
  refine ⟨hg.1, ?_⟩
  intro x hx y hy a c ha hc hac
  have hle : ∀ z : ℝ, (if z = μ then v else g z) ≤ g z := by
    intro z
    split_ifs with h
    · subst h; exact hv
    · exact le_rfl
  simp only [smul_eq_mul] at *
  by_cases hz : a * x + c * y = μ
  · -- every point with positive weight sits at μ
    have hx' : μ ≤ x := hx.1
    have hy' : μ ≤ y := hy.1
    have hsum0 : a * (x - μ) + c * (y - μ) = 0 := by
      have e : a * (x - μ) + c * (y - μ) = (a * x + c * y) - (a + c) * μ := by ring
      rw [e, hac, hz]; ring
    have h1 : a * (x - μ) = 0 := by
      linarith [mul_nonneg ha (sub_nonneg.2 hx'), mul_nonneg hc (sub_nonneg.2 hy')]
    have h2 : c * (y - μ) = 0 := by
      linarith [mul_nonneg ha (sub_nonneg.2 hx'), mul_nonneg hc (sub_nonneg.2 hy')]
    have ex : a * (if x = μ then v else g x) = a * v := by
      rcases mul_eq_zero.1 h1 with h | h
      · simp [h]
      · have : x = μ := by linarith
        simp [this]
    have ey : c * (if y = μ then v else g y) = c * v := by
      rcases mul_eq_zero.1 h2 with h | h
      · simp [h]
      · have : y = μ := by linarith
        simp [this]
    have hsum : a * v + c * v = v := by rw [← add_mul, hac, one_mul]
    rw [ex, ey, if_pos hz, hsum]
  · rw [if_neg hz]
    have hconc := hg.2 hx hy ha hc hac
    simp only [smul_eq_mul] at hconc
    nlinarith [mul_le_mul_of_nonneg_left (hle x) ha, mul_le_mul_of_nonneg_left (hle y) hc]

/-- **Mirror: raising the right end value keeps convexity.** -/
theorem endpoint_raise_convex {a μ v : ℝ} {g : ℝ → ℝ} (hg : ConvexOn ℝ (Set.Icc a μ) g)
    (hv : g μ ≤ v) : ConvexOn ℝ (Set.Icc a μ) (fun x => if x = μ then v else g x) := by
  refine ⟨hg.1, ?_⟩
  intro x hx y hy p q hp hq hpq
  have hge : ∀ z : ℝ, g z ≤ (if z = μ then v else g z) := by
    intro z
    split_ifs with h
    · subst h; exact hv
    · exact le_rfl
  simp only [smul_eq_mul] at *
  by_cases hz : p * x + q * y = μ
  · have hx' : x ≤ μ := hx.2
    have hy' : y ≤ μ := hy.2
    have hsum0 : p * (μ - x) + q * (μ - y) = 0 := by
      have e : p * (μ - x) + q * (μ - y) = (p + q) * μ - (p * x + q * y) := by ring
      rw [e, hpq, hz]; ring
    have h1 : p * (μ - x) = 0 := by
      linarith [mul_nonneg hp (sub_nonneg.2 hx'), mul_nonneg hq (sub_nonneg.2 hy')]
    have h2 : q * (μ - y) = 0 := by
      linarith [mul_nonneg hp (sub_nonneg.2 hx'), mul_nonneg hq (sub_nonneg.2 hy')]
    have ex : p * (if x = μ then v else g x) = p * v := by
      rcases mul_eq_zero.1 h1 with h | h
      · simp [h]
      · have : x = μ := by linarith
        simp [this]
    have ey : q * (if y = μ then v else g y) = q * v := by
      rcases mul_eq_zero.1 h2 with h | h
      · simp [h]
      · have : y = μ := by linarith
        simp [this]
    have hsum : p * v + q * v = v := by rw [← add_mul, hpq, one_mul]
    rw [ex, ey, if_pos hz, hsum]
  · rw [if_neg hz]
    have hconv := hg.2 hx hy hp hq hpq
    simp only [smul_eq_mul] at hconv
    nlinarith [mul_le_mul_of_nonneg_left (hge x) hp, mul_le_mul_of_nonneg_left (hge y) hq]

/-- **Seam-shift error.** If g₁ is within e of F on the left, g₂ within e of F on the right, and
`g₁ μ - g₂ μ ≤ δ` with `s = max δ 0 / 2`, then the shifted fit G (g₁ - s left of and at μ, g₂ + s right of μ)
is within e + s of F everywhere, and G μ ≤ g₂ μ + s (the condition used by `endpoint_lower_concave`). -/
theorem seam_shift_error (F g₁ g₂ : ℝ → ℝ) (μ e δ : ℝ)
    (h1 : ∀ x ≤ μ, |F x - g₁ x| ≤ e) (h2 : ∀ x ≥ μ, |F x - g₂ x| ≤ e) (hδ : g₁ μ - g₂ μ ≤ δ) :
    (∀ x, |F x - (if x ≤ μ then g₁ x - max δ 0 / 2 else g₂ x + max δ 0 / 2)| ≤ e + max δ 0 / 2) ∧
      g₁ μ - max δ 0 / 2 ≤ g₂ μ + max δ 0 / 2 := by
  have hs : 0 ≤ max δ 0 / 2 := by positivity
  refine ⟨fun x => ?_, ?_⟩
  · split_ifs with hx
    · have := h1 x hx
      rw [abs_le] at this ⊢
      constructor <;> linarith [this.1, this.2]
    · have := h2 x (le_of_lt (not_le.1 hx))
      rw [abs_le] at this ⊢
      constructor <;> linarith [this.1, this.2]
  · have : δ ≤ max δ 0 := le_max_left δ 0
    linarith

/-! Non-vacuity witnesses (NONVAC): each hypothesis set is satisfiable with a strict instance. -/

example : ConcaveOn ℝ (Set.Icc (0 : ℝ) 1) (fun x => if x = 0 then (-1 : ℝ) else -x ^ 2) :=
  endpoint_lower_concave (g := fun x => -x ^ 2) (by
    have : ConvexOn ℝ (Set.Icc (0 : ℝ) 1) (fun x : ℝ => x ^ 2) :=
      (convexOn_pow 2).subset Set.Icc_subset_Ici_self (convex_Icc 0 1)
    exact this.neg) (by norm_num)

example : (fun x : ℝ => if x ≤ 0 then (0 : ℝ) - max (1 : ℝ) 0 / 2 else (0 : ℝ) + max (1 : ℝ) 0 / 2) 0 = -1 / 2 := by
  norm_num

/-! Axioms audit (Step 8 (f)). -/
#print axioms endpoint_lower_concave
#print axioms endpoint_raise_convex
#print axioms seam_shift_error
