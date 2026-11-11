/-
Thm4Glue.lean — machine checks of the gluing and flattening steps of work/DT3-THEOREM.md §4c (Theorem 4 v0.2),
Step 8 (DT6), after the independent audit. Deterministic steps only; the stochastic limit is not formalised.

1. `glue_max_convex`  : the maximum of two convex functions on a convex set is convex (cases A, C′: G₁ := max(·,·)).
2. `glue_max_error`   : if both pieces are at most e above F, and at every point at least one of them is at least
                        F − e, then the maximum is within e of F (the error bookkeeping of the glued fit).
3. `flatten_error`    : the flattening ψ = h + χ (h m − h) with 0 ≤ χ ≤ 1 and χ = 0 away from m is within the local
                        oscillation of h (case C′ completion).
4. `flatten_const`    : ψ equals h m wherever χ = 1.
5. `tail_indicator`   : if limsup X_n ≤ T (eventually X_n ≤ T + ε for every ε > 0) and T < c, then eventually X_n ≤ c:
                        the deterministic core of reverse Fatou with {T ≥ c} (Corollary 4).
-/
import Mathlib

/-- **Max of two convex functions is convex.** -/
theorem glue_max_convex {s : Set ℝ} {f g : ℝ → ℝ} (hf : ConvexOn ℝ s f) (hg : ConvexOn ℝ s g) :
    ConvexOn ℝ s (fun x => max (f x) (g x)) :=
  hf.sup hg

/-- **Error of the glued maximum.** -/
theorem glue_max_error (F g w : ℝ → ℝ) (e : ℝ) (hg : ∀ x, g x - F x ≤ e) (hw : ∀ x, w x - F x ≤ e)
    (hlow : ∀ x, F x - e ≤ g x ∨ F x - e ≤ w x) : ∀ x, |max (g x) (w x) - F x| ≤ e := by
  intro x
  rw [abs_le]
  constructor
  · rcases hlow x with h | h
    · have : g x ≤ max (g x) (w x) := le_max_left _ _
      linarith
    · have : w x ≤ max (g x) (w x) := le_max_right _ _
      linarith
  · rcases max_cases (g x) (w x) with ⟨h, _⟩ | ⟨h, _⟩
    · rw [h]; linarith [hg x]
    · rw [h]; linarith [hw x]

/-- **Flattening error.** -/
theorem flatten_error (h χ : ℝ → ℝ) (m ρ ω : ℝ) (hρ : 0 ≤ ρ) (hosc : ∀ x, |x - m| ≤ 2 * ρ → |h m - h x| ≤ ω)
    (hχ0 : ∀ x, 0 ≤ χ x) (hχ1 : ∀ x, χ x ≤ 1) (hχout : ∀ x, 2 * ρ < |x - m| → χ x = 0) :
    ∀ x, |(h x + χ x * (h m - h x)) - h x| ≤ ω := by
  intro x
  have e : (h x + χ x * (h m - h x)) - h x = χ x * (h m - h x) := by ring
  rw [e, abs_mul, abs_of_nonneg (hχ0 x)]
  by_cases hx : |x - m| ≤ 2 * ρ
  · have h1 := hosc x hx
    have hω : 0 ≤ ω := le_trans (abs_nonneg _) h1
    calc χ x * |h m - h x| ≤ 1 * ω :=
          mul_le_mul (hχ1 x) h1 (abs_nonneg _) (by norm_num)
      _ = ω := one_mul ω
  · rw [hχout x (lt_of_not_ge hx), zero_mul]
    have := hosc m (by rw [sub_self, abs_zero]; linarith)
    exact le_trans (abs_nonneg _) this

/-- **Flattened function is constant where the cutoff is 1.** -/
theorem flatten_const (h χ : ℝ → ℝ) (m x : ℝ) (hx : χ x = 1) : h x + χ x * (h m - h x) = h m := by
  rw [hx]; ring

/-- **Deterministic core of reverse Fatou with {T ≥ c}.** -/
theorem tail_indicator (X : ℕ → ℝ) (T c : ℝ) (hlim : ∀ ε > 0, ∀ᶠ n in Filter.atTop, X n ≤ T + ε) (hT : T < c) :
    ∀ᶠ n in Filter.atTop, ¬ (c < X n) := by
  have hε : 0 < c - T := sub_pos.2 hT
  filter_upwards [hlim (c - T) hε] with n hn
  intro hc
  linarith

/-! Non-vacuity witnesses (NONVAC). -/

example : |max (1 : ℝ) 0 - 1| ≤ 1 :=
  glue_max_error (fun _ => (1 : ℝ)) (fun _ => 1) (fun _ => 0) 1 (by intro x; norm_num) (by intro x; norm_num)
    (by intro x; left; norm_num) 0

example : (fun _ : ℝ => (3 : ℝ) + 1 * (5 - 3)) 0 = 5 := by norm_num

/-! Axioms audit (Step 8 (f)). -/
#print axioms glue_max_convex
#print axioms glue_max_error
#print axioms flatten_error
#print axioms flatten_const
#print axioms tail_indicator
