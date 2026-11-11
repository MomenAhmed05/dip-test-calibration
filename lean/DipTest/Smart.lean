/-
Step 15 / DT14 WP-A: abstract (non-probabilistic) cores of Theorem T-sharp and Theorem T-sat
(`DT14-SMART.md` v0.2).

Scope: these lemmas check the order-and-limit logic only. The probabilistic inputs stay cited:
C2 (DT3 Theorem 1, liminf lower bound at fixed p), continuity/strict monotonicity of the
Balabdaoui-Pitman law, the diagonal construction, Markov's inequality and contiguity transfer.

* `sharp_core`: a threshold test whose rejection probability is antitone in the threshold,
  eventually at most α + ε, while the probability at a fixed threshold b is eventually ≥ q > α,
  must eventually use thresholds above b.
* `sharp_liminf`: if such a law exists for every b ∈ [0, c), the thresholds eventually exceed
  every b < c (i.e. liminf a_n ≥ c).
* `sat_t`, `sat_core`: the numerical saturation core of T-sat.
-/
import Mathlib

open Filter Topology

namespace DipSmart

/-- T-sharp core at one threshold `b`. `P n x` = rejection probability at threshold `x`
under one fixed law (antitone in `x`). -/
theorem sharp_core (a : ℕ → ℝ) (P : ℕ → ℝ → ℝ) (α b q : ℝ)
    (hanti : ∀ n x y, x ≤ y → P n y ≤ P n x)
    (hvalid : ∀ ε > 0, ∀ᶠ n in atTop, P n (a n) ≤ α + ε)
    (hlow : ∀ᶠ n in atTop, q ≤ P n b) (hq : α < q) :
    ∀ᶠ n in atTop, b < a n := by
  have hε : (0 : ℝ) < (q - α) / 2 := by linarith
  filter_upwards [hvalid _ hε, hlow] with n h1 h2
  by_contra hcon
  have h3 := hanti n (a n) b (not_lt.1 hcon)
  linarith

/-- T-sharp core, all thresholds below `c`: for every `b ∈ [0, c)` some law (its own
`P`, `q`) witnesses the hypotheses of `sharp_core`; thresholds are nonnegative. -/
theorem sharp_liminf (a : ℕ → ℝ) (α c : ℝ) (ha0 : ∀ n, 0 ≤ a n)
    (hfam : ∀ b, 0 ≤ b → b < c → ∃ (P : ℕ → ℝ → ℝ) (q : ℝ),
      (∀ n x y, x ≤ y → P n y ≤ P n x) ∧
      (∀ ε > 0, ∀ᶠ n in atTop, P n (a n) ≤ α + ε) ∧
      (∀ᶠ n in atTop, q ≤ P n b) ∧ α < q) :
    ∀ b < c, ∀ᶠ n in atTop, b < a n := by
  intro b hb
  rcases lt_or_ge b 0 with hneg | hnn
  · exact Eventually.of_forall (fun n => lt_of_lt_of_le hneg (ha0 n))
  · obtain ⟨P, q, hanti, hvalid, hlow, hq⟩ := hfam b hnn hb
    exact sharp_core a P α b q hanti hvalid hlow hq

/-- T-sat: the null rejection probability of the fix along the saturating sequence tends to α. -/
theorem sat_t (s t : ℕ → ℝ) (α : ℝ) (hts : ∀ n, t n ≤ s n)
    (hs : ∀ ε > 0, ∀ᶠ n in atTop, s n ≤ α + ε)
    (ht : ∀ ε > 0, ∀ᶠ n in atTop, α - ε ≤ t n) :
    Tendsto t atTop (𝓝 α) := by
  rw [Metric.tendsto_atTop]
  intro ε hε
  have h2 : (0 : ℝ) < ε / 2 := by positivity
  obtain ⟨N1, hN1⟩ := eventually_atTop.1 (hs (ε / 2) h2)
  obtain ⟨N2, hN2⟩ := eventually_atTop.1 (ht (ε / 2) h2)
  refine ⟨max N1 N2, fun n hn => ?_⟩
  have a1 := hN1 n (le_of_max_le_left hn)
  have a2 := hN2 n (le_of_max_le_right hn)
  have a3 := hts n
  rw [Real.dist_eq, abs_lt]
  constructor <;> linarith

/-- T-sat numerical core: 0 ≤ e_n ≤ s_n − t_n with limsup s_n ≤ α and liminf t_n ≥ α
(in eventually-ε form) forces e_n → 0. -/
theorem sat_core (e s t : ℕ → ℝ) (α : ℝ)
    (h0 : ∀ n, 0 ≤ e n) (h1 : ∀ n, e n ≤ s n - t n)
    (hs : ∀ ε > 0, ∀ᶠ n in atTop, s n ≤ α + ε)
    (ht : ∀ ε > 0, ∀ᶠ n in atTop, α - ε ≤ t n) :
    Tendsto e atTop (𝓝 0) := by
  rw [Metric.tendsto_atTop]
  intro ε hε
  have h2 : (0 : ℝ) < ε / 4 := by positivity
  obtain ⟨N1, hN1⟩ := eventually_atTop.1 (hs (ε / 4) h2)
  obtain ⟨N2, hN2⟩ := eventually_atTop.1 (ht (ε / 4) h2)
  refine ⟨max N1 N2, fun n hn => ?_⟩
  have a1 := hN1 n (le_of_max_le_left hn)
  have a2 := hN2 n (le_of_max_le_right hn)
  have a3 := h1 n
  rw [Real.dist_eq, sub_zero, abs_of_nonneg (h0 n)]
  linarith

/-- Contiguity-transfer bound used in T-sat (ii): for Δ ≤ 1 and ε ≥ 0,
E Δ ≤ ε + Q(Δ > ε), stated for a finite probability vector (a non-vacuous algebraic core). -/
theorem expect_le_eps_add_tail {ι : Type*} (S : Finset ι) (w Δ : ι → ℝ) (ε : ℝ) (hε : 0 ≤ ε)
    (hw : ∀ i ∈ S, 0 ≤ w i) (hsum : ∑ i ∈ S, w i = 1)
    (hΔ1 : ∀ i ∈ S, Δ i ≤ 1) :
    ∑ i ∈ S, w i * Δ i ≤ ε + ∑ i ∈ S.filter (fun i => ε < Δ i), w i := by
  classical
  have key : ∀ i ∈ S, w i * Δ i ≤ w i * ε + (if ε < Δ i then w i else 0) := by
    intro i hi
    by_cases h : ε < Δ i
    · simp only [h, ite_true]
      have h4 : w i * Δ i ≤ w i * 1 := mul_le_mul_of_nonneg_left (hΔ1 i hi) (hw i hi)
      have h5 : 0 ≤ w i * ε := mul_nonneg (hw i hi) hε
      linarith
    · simp only [h, ite_false, add_zero]
      exact mul_le_mul_of_nonneg_left (not_lt.1 h) (hw i hi)
  calc ∑ i ∈ S, w i * Δ i ≤ ∑ i ∈ S, (w i * ε + (if ε < Δ i then w i else 0)) :=
        Finset.sum_le_sum key
    _ = ε + ∑ i ∈ S.filter (fun i => ε < Δ i), w i := by
        rw [Finset.sum_add_distrib, ← Finset.sum_mul, hsum, one_mul, Finset.sum_filter]

/-- Non-vacuity witness for `sharp_core`. -/
example : ∀ᶠ n in atTop, (0 : ℝ) < (fun _ : ℕ => (1 : ℝ)) n :=
  sharp_core (fun _ => 1) (fun _ x => if x < 1 then 1 else 0) (1 / 2) 0 1
    (by intro n x y hxy; by_cases hy : y < 1 <;> by_cases hx : x < 1 <;> simp [hx, hy] <;> linarith)
    (by intro ε hε; exact Eventually.of_forall (fun n => by simp; linarith))
    (Eventually.of_forall (fun n => by simp)) (by norm_num)

end DipSmart

#print axioms DipSmart.sharp_core
#print axioms DipSmart.sharp_liminf
#print axioms DipSmart.sat_t
#print axioms DipSmart.sat_core
#print axioms DipSmart.expect_le_eps_add_tail
