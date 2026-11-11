/-
DT11 Step 12 (work/DT11-NEAREDGE.md): machine checks for the near-edge band of C-edge∞.

Main results, for t ∈ [0.72, 0.88] and ε ∈ (0, 0.15]:
* `near_edge_band_published` (Theorem NE end to end): g < 1 − F(2t) from the published inputs
  P1–P5 for F = F_M and the Step-11 Lemma B for g.
* `lemma1_increments` (Lemma 1, increment form): P1, P2, P5 ⇒ UP: F d − F c ≤ ∫_c^d I₁ on
  [1, 1.8] and LO: 0.935·∫_c^d I₁ ≤ F d − F c on [1.44, 1.8].
* `near_edge_band` (Theorem NE): UP, LO + the bridge inputs + Lemma B ⇒ g < 1 − F(2t);
  `near_edge_band_density` is the same theorem with a pointwise density sandwich.

What Lean checks here (our new near-edge reasoning): Lemma 1 from its published inputs (finite
product rule and FTC per stick sequence, Weierstrass product inequality, Tonelli, the change of
variables v = z/√ℓ, the tail integral ∫_z^∞ f_E = S(z), the first-moment identity E Σ h_j = I₁,
the ordered-pair bound with J ≤ J̄, monotone limits in n); the series bounds for S = 1 − F_E
(Lemma 2); the uniform box inequality (`scalar_box`); the layer-cake gain bound
(`expect_exp_min_le`, `gain_bound`); the loss bound; the assembly; the path facts of Lemma 3
(`gcm_gap_le_range`, `Lambda_le_X`); and a non-vacuity witness for `near_edge_band`
(`near_edge_band_nonvac`).

What enters as HYPOTHESES (published results are not re-proved):
  P1 Balabdaoui–Pitman 2011 + independence of L and (E_j) + {M_n ≤ z} ↓ {M ≤ z}
     (`hP1`: F z = lim_n E Π_{i<n} F_E(z/√L_i); F is a CDF),
  P2 Kennedy 1976: F_E = P is a CDF, P = 1 − S on [1, ∞), with a continuous density there,
  P3 Kuiper/Vervaat range law (μ(Aᶜ), μ(E) ≤ S(0.6t/ε)),
  P4 bridge linear-boundary crossing (`hcross`),
  P5 Pitman: GEM(1) first and ordered-pair moment measures (`hP5a`, `hP5b`),
  and the Step-11 Lemma B (`hB`). Identifying F, g, T, X, L with the Brownian/GEM objects is the
  paper's instantiation step: Mathlib has no Brownian-bridge path theory and no
  GEM(1) law. The witness shows the hypotheses of `near_edge_band` are jointly satisfiable; it does
  not show that the BP11 law satisfies them (that is P1–P5, published).
No `sorry`, no custom axioms (every declaration: [propext, Classical.choice, Quot.sound]).
Compile: `cd ../rough-vol/lean && lake env lean <dip-test>/lean/NearEdge.lean`.
-/
import Mathlib.MeasureTheory.Integral.Layercake
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic
import Mathlib.Analysis.Complex.ExponentialBounds
import Mathlib.Analysis.SpecificLimits.Basic
import Mathlib.Analysis.Convex.Function
import Mathlib.Analysis.SpecialFunctions.Sqrt
import Mathlib.Analysis.SpecialFunctions.ImproperIntegrals
import Mathlib.MeasureTheory.Function.JacobianOneDim
import Mathlib.MeasureTheory.Integral.IntegralEqImproper
import Mathlib.Analysis.Calculus.Deriv.Slope

open MeasureTheory Set Real

open scoped ENNReal

namespace NearEdge

noncomputable def term (z : ℝ) (k : ℕ) : ℝ :=
  (4 * ((k : ℝ) + 1) ^ 2 * z ^ 2 - 1) * exp (-2 * ((k : ℝ) + 1) ^ 2 * z ^ 2)

noncomputable def S (z : ℝ) : ℝ := 2 * ∑' k : ℕ, term z k

noncomputable def Jbar (z : ℝ) : ℝ := 4 * ∑' k : ℕ, exp (-2 * ((k : ℝ) + 1) ^ 2 * z ^ 2)

/-- e > 2.7182818283 gives e^n lower bounds. -/
lemma exp_nat_ge (n : ℕ) : (2.7182818283 : ℝ) ^ n ≤ exp n := by
  rw [← Real.exp_one_pow]
  exact pow_le_pow_left₀ (by norm_num) Real.exp_one_gt_d9.le n

lemma exp_nat_le (n : ℕ) : exp n ≤ (2.7182818286 : ℝ) ^ n := by
  rw [← Real.exp_one_pow]
  exact pow_le_pow_left₀ (Real.exp_pos 1).le Real.exp_one_lt_d9.le n

lemma exp_neg_lt_one {c : ℝ} (hc : 0 < c) : exp (-c) < 1 := by
  have := Real.exp_lt_exp.mpr (show -c < 0 by linarith)
  rwa [Real.exp_zero] at this

lemma term_nonneg {z : ℝ} (hz : 1 ≤ z) (k : ℕ) : 0 ≤ term z k := by
  unfold term
  have hk : (1 : ℝ) ≤ (k : ℝ) + 1 := by have := Nat.cast_nonneg (α := ℝ) k; linarith
  have h1 : (1 : ℝ) ≤ ((k : ℝ) + 1) ^ 2 := one_le_pow₀ hk
  have h2 : (1 : ℝ) ≤ z ^ 2 := one_le_pow₀ hz
  have : 0 ≤ 4 * ((k : ℝ) + 1) ^ 2 * z ^ 2 - 1 := by nlinarith
  positivity

/-- Geometric domination of the higher terms: for z ≥ 1,
term z (k+1) ≤ term z 0 · (16/3) e^{-6} (e^{-5})^k. -/
lemma term_succ_le {z : ℝ} (hz : 1 ≤ z) (k : ℕ) :
    term z (k + 1) ≤ term z 0 * ((16 / 3) * exp (-6) * (exp (-5)) ^ k) := by
  unfold term
  push_cast
  set n : ℝ := (k : ℝ) + 2 with hn
  have hk0 : (0 : ℝ) ≤ k := Nat.cast_nonneg k
  have hz2 : (1 : ℝ) ≤ z ^ 2 := one_le_pow₀ hz
  -- polynomial factor: 4n²z² − 1 ≤ (4/3) n² (4z² − 1)
  have hpoly : 4 * ((k : ℝ) + 1 + 1) ^ 2 * z ^ 2 - 1 ≤ (4 / 3) * n ^ 2 * (4 * z ^ 2 - 1) := by
    rw [hn]; nlinarith [sq_nonneg ((k : ℝ) + 2)]
  -- exponential factor
  have hexp : exp (-2 * ((k : ℝ) + 1 + 1) ^ 2 * z ^ 2)
      ≤ exp (-2 * z ^ 2) * exp (-2 * (n ^ 2 - 1)) := by
    rw [← Real.exp_add]; apply Real.exp_le_exp.mpr; rw [hn]
    have hn1 : 0 ≤ ((k : ℝ) + 2) ^ 2 - 1 := by nlinarith
    nlinarith [mul_nonneg hn1 (show 0 ≤ z ^ 2 - 1 by linarith)]
  -- n² e^{-2(n²-1)} ≤ 4 e^{-6} e^{-5k}
  have hgeo : n ^ 2 * exp (-2 * (n ^ 2 - 1)) ≤ 4 * exp (-6) * (exp (-5)) ^ k := by
    have hy : n ^ 2 / 4 ≤ exp ((2 * (k : ℝ) + 3) * k) := by
      have := Real.add_one_le_exp ((2 * (k : ℝ) + 3) * k)
      rw [hn]; nlinarith
    rw [← Real.exp_nat_mul]
    have he : exp (-2 * (n ^ 2 - 1)) * exp ((2 * (k : ℝ) + 3) * k) = exp (-6) * exp (k * -5) := by
      rw [← Real.exp_add, ← Real.exp_add]; congr 1; rw [hn]; ring
    have hpos : 0 < exp (-2 * (n ^ 2 - 1)) := Real.exp_pos _
    calc n ^ 2 * exp (-2 * (n ^ 2 - 1)) = 4 * (n ^ 2 / 4) * exp (-2 * (n ^ 2 - 1)) := by ring
      _ ≤ 4 * exp ((2 * (k : ℝ) + 3) * k) * exp (-2 * (n ^ 2 - 1)) := by gcongr
      _ = 4 * (exp (-2 * (n ^ 2 - 1)) * exp ((2 * (k : ℝ) + 3) * k)) := by ring
      _ = 4 * exp (-6) * exp (k * -5) := by rw [he]; ring
  have hA0 : 0 ≤ 4 * z ^ 2 - 1 := by nlinarith
  have hE0 : 0 ≤ exp (-2 * z ^ 2) := (Real.exp_pos _).le
  have hlhs0 : 0 ≤ 4 * ((k : ℝ) + 1 + 1) ^ 2 * z ^ 2 - 1 := by nlinarith [sq_nonneg ((k : ℝ) + 2)]
  calc (4 * ((k : ℝ) + 1 + 1) ^ 2 * z ^ 2 - 1) * exp (-2 * ((k : ℝ) + 1 + 1) ^ 2 * z ^ 2)
      ≤ ((4 / 3) * n ^ 2 * (4 * z ^ 2 - 1)) * (exp (-2 * z ^ 2) * exp (-2 * (n ^ 2 - 1))) := by
        apply mul_le_mul hpoly hexp (Real.exp_pos _).le; positivity
    _ = (4 * z ^ 2 - 1) * exp (-2 * z ^ 2) * ((4 / 3) * (n ^ 2 * exp (-2 * (n ^ 2 - 1)))) := by ring
    _ ≤ (4 * z ^ 2 - 1) * exp (-2 * z ^ 2) * ((4 / 3) * (4 * exp (-6) * (exp (-5)) ^ k)) := by
        gcongr
    _ = (4 * (0 + 1) ^ 2 * z ^ 2 - 1) * exp (-2 * (0 + 1) ^ 2 * z ^ 2)
          * ((16 / 3) * exp (-6) * (exp (-5)) ^ k) := by ring

lemma summable_term {z : ℝ} (hz : 1 ≤ z) : Summable (term z) := by
  have hg : Summable (fun k : ℕ => term z 0 * ((16 / 3) * exp (-6) * (exp (-5)) ^ k)) := by
    apply Summable.mul_left; apply Summable.mul_left
    exact summable_geometric_of_lt_one (Real.exp_pos _).le (exp_neg_lt_one (by norm_num))
  have h1 : Summable (fun k => term z (k + 1)) :=
    Summable.of_nonneg_of_le (fun k => term_nonneg hz (k + 1)) (term_succ_le hz) hg
  exact (summable_nat_add_iff 1).mp h1

theorem S_lower {z : ℝ} (hz : 1 ≤ z) : 2 * (4 * z ^ 2 - 1) * exp (-2 * z ^ 2) ≤ S z := by
  unfold S
  have h := (summable_term hz).le_tsum 0 (fun j _ => term_nonneg hz j)
  have h0 : term z 0 = (4 * z ^ 2 - 1) * exp (-2 * z ^ 2) := by unfold term; push_cast; ring_nf
  rw [h0] at h; linarith

theorem S_upper {z : ℝ} (hz : 1 ≤ z) :
    S z ≤ 2 * (4 * z ^ 2 - 1) * exp (-2 * z ^ 2) * (1 + 0.014) := by
  unfold S
  have hs := summable_term hz
  rw [hs.tsum_eq_zero_add]
  have hg : Summable (fun k : ℕ => term z 0 * ((16 / 3) * exp (-6) * (exp (-5)) ^ k)) := by
    apply Summable.mul_left; apply Summable.mul_left
    exact summable_geometric_of_lt_one (Real.exp_pos _).le (exp_neg_lt_one (by norm_num))
  have htail : ∑' k, term z (k + 1) ≤ ∑' k : ℕ, term z 0 * ((16 / 3) * exp (-6) * (exp (-5)) ^ k) :=
    Summable.tsum_le_tsum (term_succ_le hz) ((summable_nat_add_iff 1).mpr hs) hg
  have hsum : ∑' k : ℕ, term z 0 * ((16 / 3) * exp (-6) * (exp (-5)) ^ k)
      = term z 0 * ((16 / 3) * exp (-6) * (1 - exp (-5))⁻¹) := by
    rw [tsum_mul_left, tsum_mul_left, tsum_geometric_of_lt_one (Real.exp_pos _).le
      (exp_neg_lt_one (by norm_num))]
  -- numeric: (16/3) e^{-6} / (1 - e^{-5}) ≤ 0.014
  have he6 : exp (-6) ≤ 1 / 403 := by
    rw [Real.exp_neg, inv_eq_one_div]
    apply one_div_le_one_div_of_le (by norm_num)
    have := exp_nat_ge 6; push_cast at this; norm_num at this ⊢; linarith
  have he5 : exp (-5) ≤ 1 / 148 := by
    rw [Real.exp_neg, inv_eq_one_div]
    apply one_div_le_one_div_of_le (by norm_num)
    have := exp_nat_ge 5; push_cast at this; norm_num at this ⊢; linarith
  have hc : (16 / 3) * exp (-6) * (1 - exp (-5))⁻¹ ≤ 0.014 := by
    have h1 : 0 < 1 - exp (-5) := by linarith
    rw [← div_eq_mul_inv, div_le_iff₀ h1]; nlinarith
  have h0 : term z 0 = (4 * z ^ 2 - 1) * exp (-2 * z ^ 2) := by unfold term; push_cast; ring_nf
  have ht0 : 0 ≤ term z 0 := term_nonneg hz 0
  rw [hsum] at htail
  have : term z 0 * ((16 / 3) * exp (-6) * (1 - exp (-5))⁻¹) ≤ term z 0 * 0.014 :=
    mul_le_mul_of_nonneg_left hc ht0
  rw [h0] at htail this ⊢
  linarith

theorem Jbar_le {z : ℝ} (hz : 1.44 ≤ z) : Jbar z ≤ 0.065 := by
  unfold Jbar
  set q : ℝ := exp (-4.1472) with hq
  have hq0 : 0 ≤ q := (Real.exp_pos _).le
  have hq1 : q < 1 := by rw [hq]; exact exp_neg_lt_one (by norm_num)
  have hz2 : (2.0736 : ℝ) ≤ z ^ 2 := by nlinarith
  have hle : ∀ k : ℕ, exp (-2 * ((k : ℝ) + 1) ^ 2 * z ^ 2) ≤ q * q ^ k := by
    intro k
    rw [hq, ← Real.exp_nat_mul, ← Real.exp_add]
    apply Real.exp_le_exp.mpr
    have hk0 : (0 : ℝ) ≤ k := Nat.cast_nonneg k
    have : ((k : ℝ) + 1) ≤ ((k : ℝ) + 1) ^ 2 := by nlinarith
    nlinarith
  have hg : Summable (fun k : ℕ => q * q ^ k) :=
    (summable_geometric_of_lt_one hq0 hq1).mul_left q
  have hs : Summable (fun k : ℕ => exp (-2 * ((k : ℝ) + 1) ^ 2 * z ^ 2)) :=
    Summable.of_nonneg_of_le (fun k => (Real.exp_pos _).le) hle hg
  have ht := Summable.tsum_le_tsum hle hs hg
  rw [tsum_mul_left, tsum_geometric_of_lt_one hq0 hq1] at ht
  -- q ≤ 1/63
  have hqv : q ≤ 1 / 63 := by
    rw [hq, show (-4.1472 : ℝ) = -(4 + 0.1472) by norm_num, Real.exp_neg, inv_eq_one_div]
    apply one_div_le_one_div_of_le (by norm_num)
    rw [Real.exp_add]
    have h4 := exp_nat_ge 4; push_cast at h4
    have h01 := Real.quadratic_le_exp_of_nonneg (x := 0.1472) (by norm_num)
    have : (63 : ℝ) ≤ 2.7182818283 ^ 4 * (1 + 0.1472 + 0.1472 ^ 2 / 2) := by norm_num
    nlinarith [Real.exp_pos 4, Real.exp_pos 0.1472]
  have h1q : 0 < 1 - q := by linarith
  have : q * (1 - q)⁻¹ ≤ 1 / 62 := by
    rw [← div_eq_mul_inv, div_le_iff₀ h1q]; nlinarith
  nlinarith


/-- Facts about a = √(1-ε²) on the box. -/
lemma box_a (ε a : ℝ) (hε0 : 0 < ε) (hε1 : ε ≤ 0.15) (ha : a = √(1 - ε ^ 2)) :
    a ^ 2 = 1 - ε ^ 2 ∧ 0 < a ∧ a ≤ 1 ∧ 0.9886 ≤ a ∧ ε ^ 2 ≤ 0.0225 := by
  have he2 : ε ^ 2 ≤ 0.0225 := by nlinarith
  have ha2 : a ^ 2 = 1 - ε ^ 2 := by rw [ha, Real.sq_sqrt]; nlinarith
  have ha0 : 0 < a := by rw [ha]; apply Real.sqrt_pos.mpr; nlinarith
  exact ⟨ha2, ha0, by nlinarith, by nlinarith, he2⟩

lemma box_E2 (x : ℝ) (hx2 : x ^ 2 ≤ 3.17) : 0.0016 ≤ exp (-2 * x ^ 2) := by
  have h1 : exp (-6.34) ≤ exp (-2 * x ^ 2) := Real.exp_le_exp.mpr (by nlinarith)
  have h6 : exp 6 ≤ 2.7182818286 ^ 6 := by have := exp_nat_le 6; push_cast at this; exact this
  have h34 : 1 - 0.34 ≤ exp (-0.34) := by have := Real.add_one_le_exp (-0.34); linarith
  have hprod : exp (-6.34) * (exp 6 * (1 / exp (-0.34))) = 1 := by
    rw [show exp (-6.34) * (exp 6 * (1 / exp (-0.34))) = exp (-6.34 + 6) / exp (-0.34) by
      rw [Real.exp_add]; ring]
    rw [show (-6.34 : ℝ) + 6 = -0.34 by norm_num, div_self (Real.exp_pos _).ne']
  have hinv : 1 / exp (-0.34) ≤ 1 / 0.66 := one_div_le_one_div_of_le (by norm_num) (by linarith)
  have hbig : exp 6 * (1 / exp (-0.34)) ≤ 625 := by
    have : (2.7182818286 : ℝ) ^ 6 * (1 / 0.66) ≤ 625 := by norm_num
    calc exp 6 * (1 / exp (-0.34)) ≤ 2.7182818286 ^ 6 * (1 / 0.66) := by
          apply mul_le_mul h6 hinv (by positivity) (by positivity)
      _ ≤ 625 := this
  have hpos : 0 < exp 6 * (1 / exp (-0.34)) := by positivity
  nlinarith [Real.exp_pos (-6.34)]

lemma box_loss (t ε a x E2 : ℝ) (ht1 : 0.72 ≤ t) (hε0 : 0 < ε) (ha2 : a ^ 2 = 1 - ε ^ 2)
    (ha0 : 0 < a) (hxa : x * a = 2 * t) (hE2 : 0 < E2) :
    13.6 * ε ^ 2 * E2 ≤ (x - 2 * t) * (4 * (8 * t - 1 / (2 * t)) * E2 * (1 - 0.065)) := by
  have ht0 : 0 < t := by linarith
  have he2p : 0 < ε ^ 2 := by positivity
  have hag : a * (2 + ε ^ 2) ≤ 2 := by
    nlinarith [sq_nonneg (ε ^ 2), mul_pos ha0 he2p]
  have hxt : t * ε ^ 2 ≤ x - 2 * t := by
    have : (x - 2 * t) * a = 2 * t - 2 * t * a := by rw [sub_mul, hxa]
    nlinarith
  have hm : 4 * (8 * t - 1 / (2 * t)) * E2 * (1 - 0.065) = 3.74 * (8 * t ^ 2 - 0.5) / t * E2 := by
    field_simp; ring
  rw [hm]
  have hc : 13.6 ≤ 3.74 * (8 * t ^ 2 - 0.5) := by nlinarith
  have hpos : 0 ≤ 3.74 * (8 * t ^ 2 - 0.5) / t * E2 := by
    apply mul_nonneg _ hE2.le; apply div_nonneg _ ht0.le; linarith
  have heq : t * ε ^ 2 * (3.74 * (8 * t ^ 2 - 0.5) / t * E2)
      = ε ^ 2 * E2 * (3.74 * (8 * t ^ 2 - 0.5)) := by
    field_simp
  have h1 : 13.6 * ε ^ 2 * E2 ≤ t * ε ^ 2 * (3.74 * (8 * t ^ 2 - 0.5) / t * E2) := by
    rw [heq]
    have := mul_le_mul_of_nonneg_left hc (mul_pos he2p hE2).le
    linarith
  exact h1.trans (mul_le_mul_of_nonneg_right hxt hpos)

lemma box_gain (t ε a x δ κ E2 : ℝ) (ht1 : 0.72 ≤ t) (ht2 : t ≤ 0.88) (hε0 : 0 < ε)
    (ha2 : a ^ 2 = 1 - ε ^ 2) (ha0 : 0 < a) (he2 : ε ^ 2 ≤ 0.0225)
    (hx : x = 2 * t / a) (hδ : δ = ε / a) (hκ : κ = 1.4 * t / ε) (hE2 : 0 < E2) :
    0 < 2 * κ - 4 * x * δ ∧
    δ * (4 * (4 * x - 1 / x) * (1 + 0.014) * E2) / (2 * κ - 4 * x * δ) ≤ 11.7 * ε ^ 2 * E2 := by
  have ht0 : 0 < t := by linarith
  have he2p : 0 < ε ^ 2 := by positivity
  have hD : 2 * κ - 4 * x * δ = t * (2.8 * a ^ 2 - 8 * ε ^ 2) / (ε * a ^ 2) := by
    rw [hκ, hδ, hx]; field_simp; ring
  have hDp : 0 < 2.8 * a ^ 2 - 8 * ε ^ 2 := by nlinarith
  have hDpos : 0 < 2 * κ - 4 * x * δ := by rw [hD]; positivity
  refine ⟨hDpos, ?_⟩
  have hN : δ * (4 * (4 * x - 1 / x) * (1 + 0.014) * E2)
      = ε * E2 * (4.056 * (16 * t ^ 2 - a ^ 2)) / (2 * t * a ^ 2) := by
    rw [hδ, hx]; field_simp; ring
  have hcore : 4.056 * (16 * t ^ 2 - a ^ 2) ≤ 23.4 * t ^ 2 * (2.8 * a ^ 2 - 8 * ε ^ 2) := by
    rw [ha2]
    have hT1 : 0.5184 ≤ t ^ 2 := by nlinarith
    have hT2 : t ^ 2 ≤ 0.7744 := by nlinarith
    nlinarith [mul_nonneg (sub_nonneg.mpr hT2) (sub_nonneg.mpr he2),
      mul_nonneg (sub_nonneg.mpr hT1) he2p.le]
  rw [hN, div_le_iff₀ hDpos, hD]
  have hlhs : ε * E2 * (4.056 * (16 * t ^ 2 - a ^ 2)) / (2 * t * a ^ 2)
      = (ε * E2 / (2 * t * a ^ 2)) * (4.056 * (16 * t ^ 2 - a ^ 2)) := by ring
  have hrhs : 11.7 * ε ^ 2 * E2 * (t * (2.8 * a ^ 2 - 8 * ε ^ 2) / (ε * a ^ 2))
      = (ε * E2 / (2 * t * a ^ 2)) * (23.4 * t ^ 2 * (2.8 * a ^ 2 - 8 * ε ^ 2)) := by
    field_simp; ring
  rw [hlhs, hrhs]
  exact mul_le_mul_of_nonneg_left hcore (by positivity)

lemma box_tail1 (t ε : ℝ) (ht1 : 0.72 ≤ t) (hε0 : 0 < ε) (hε1 : ε ≤ 0.15) :
    2 * S (0.6 * t / ε) ≤ 0.001 * ε ^ 2 := by
  set w := 0.6 * t / ε with hw
  have hwe : w * ε = 0.6 * t := by rw [hw]; field_simp
  have hw1 : 2.88 ≤ w := by
    rw [hw, le_div_iff₀ hε0]; nlinarith
  have hS := S_upper (show (1 : ℝ) ≤ w by linarith)
  set s := w ^ 2 with hs
  have hs8 : 8 ≤ s := by nlinarith
  have hse : s * ε ^ 2 = 0.36 * t ^ 2 := by rw [hs, ← mul_pow, hwe]; ring
  have hsq : s ^ 2 * exp (-2 * s) ≤ 64 * exp (-16) := by
    have h1 := Real.add_one_le_exp (s / 8 - 1)
    have h2 : s ≤ 8 * exp (s / 8 - 1) := by linarith
    have h3 : s ^ 2 ≤ 64 * exp (s / 4 - 2) := by
      have : exp (s / 8 - 1) ^ 2 = exp (s / 4 - 2) := by rw [← Real.exp_nat_mul]; ring_nf
      nlinarith [Real.exp_pos (s / 8 - 1)]
    have h4 : exp (s / 4 - 2) * exp (-2 * s) ≤ exp (-16) := by
      rw [← Real.exp_add]; apply Real.exp_le_exp.mpr; linarith
    nlinarith [Real.exp_pos (-2 * s), Real.exp_pos (s / 4 - 2)]
  have he16 : exp (-16) ≤ 1 / 8800000 := by
    rw [Real.exp_neg, inv_eq_one_div]
    apply one_div_le_one_div_of_le (by norm_num)
    have := exp_nat_ge 16; push_cast at this; norm_num at this ⊢; linarith
  have hSw : 2 * S w ≤ 16.224 * s * exp (-2 * s) := by
    have : exp (-2 * w ^ 2) = exp (-2 * s) := by rw [hs]
    rw [this] at hS; nlinarith [Real.exp_pos (-2 * s)]
  have hkey : 16.224 * s ^ 2 * exp (-2 * s) ≤ 0.001 * (0.36 * t ^ 2) := by nlinarith
  have hspos : 0 < s := by linarith
  have h5 : 16.224 * s * exp (-2 * s) ≤ 0.001 * ε ^ 2 := by
    have h := hkey; rw [← hse] at h
    have h6 : s * (16.224 * s * exp (-2 * s)) ≤ s * (0.001 * ε ^ 2) := by nlinarith
    exact le_of_mul_le_mul_left h6 hspos
  linarith

lemma box_tail2 (t ε a R κ : ℝ) (ht1 : 0.72 ≤ t) (hε0 : 0 < ε) (he2 : ε ^ 2 ≤ 0.0225)
    (ha1 : a ≤ 1) (hRe : R * ε = 2 * t - a) (hκe : κ * ε = 1.4 * t) :
    0 < R ∧ exp (-2 * R * (R + κ)) ≤ 0.0001 * ε ^ 2 := by
  have he2p : 0 < ε ^ 2 := by positivity
  have hR0 : 0.44 ≤ R * ε := by rw [hRe]; linarith
  have hRpos : 0 < R := by
    by_contra! h; nlinarith
  refine ⟨hRpos, ?_⟩
  set σ := 1 / ε ^ 2 with hσ
  have hσe : σ * ε ^ 2 = 1 := by rw [hσ]; field_simp
  have hσ44 : 44 ≤ σ := by rw [hσ, le_div_iff₀ he2p]; nlinarith
  have hexp : -2 * R * (R + κ) ≤ -1.27 * σ := by
    have h2 : 2 * R * (R + κ) * ε ^ 2 ≥ 1.27 := by
      have : 2 * R * (R + κ) * ε ^ 2 = 2 * (R * ε) * (R * ε + κ * ε) := by ring
      rw [this, hκe]; nlinarith
    nlinarith
  have hσexp : σ ≤ exp σ := by have := Real.add_one_le_exp σ; linarith
  have h1 : σ * exp (-1.27 * σ) ≤ exp (-0.27 * σ) := by
    calc σ * exp (-1.27 * σ) ≤ exp σ * exp (-1.27 * σ) := by gcongr
      _ = exp (-0.27 * σ) := by rw [← Real.exp_add]; ring_nf
  have h2 : exp (-0.27 * σ) ≤ exp (-11) := Real.exp_le_exp.mpr (by nlinarith)
  have h3 : exp (-11) ≤ 0.0001 := by
    rw [Real.exp_neg, inv_eq_one_div, div_le_iff₀ (Real.exp_pos _)]
    have := exp_nat_ge 11; push_cast at this; nlinarith
  have hσpos : 0 < σ := by linarith
  have h4 : σ * exp (-2 * R * (R + κ)) ≤ 0.0001 := by
    have := Real.exp_le_exp.mpr hexp
    nlinarith
  have h5 : exp (-2 * R * (R + κ)) = ε ^ 2 * (σ * exp (-2 * R * (R + κ))) := by
    rw [← mul_assoc, mul_comm (ε ^ 2), hσe, one_mul]
  rw [h5]; nlinarith

/-- The uniform box inequality (DT11 §5 steps 5-6). -/
lemma scalar_box (t ε a x δ κ R : ℝ) (ht1 : 0.72 ≤ t) (ht2 : t ≤ 0.88) (hε0 : 0 < ε)
    (hε1 : ε ≤ 0.15) (ha : a = √(1 - ε ^ 2)) (hx : x = 2 * t / a) (hδ : δ = ε / a)
    (hκ : κ = 1.4 * t / ε) (hR : R = (x - 1) / δ) :
    0 < R ∧ 0 < 2 * κ - 4 * x * δ ∧
    2 * S (0.6 * t / ε)
      + δ * (4 * (4 * x - 1 / x) * (1 + 0.014) * exp (-2 * x ^ 2)) / (2 * κ - 4 * x * δ)
      + exp (-2 * R * (R + κ))
      < (x - 2 * t) * (4 * (8 * t - 1 / (2 * t)) * exp (-2 * x ^ 2) * (1 - 0.065)) := by
  obtain ⟨ha2, ha0, ha1, _haL, he2⟩ := box_a ε a hε0 hε1 ha
  have hxa : x * a = 2 * t := by rw [hx]; field_simp
  have hx2 : x ^ 2 ≤ 3.17 := by
    have : x ^ 2 * a ^ 2 = 4 * t ^ 2 := by rw [← mul_pow, hxa]; ring
    nlinarith
  have hE2 := box_E2 x hx2
  have hE2pos : 0 < exp (-2 * x ^ 2) := Real.exp_pos _
  have hl := box_loss t ε a x (exp (-2 * x ^ 2)) ht1 hε0 ha2 ha0 hxa hE2pos
  obtain ⟨hD, hg⟩ :=
    box_gain t ε a x δ κ (exp (-2 * x ^ 2)) ht1 ht2 hε0 ha2 ha0 he2 hx hδ hκ hE2pos
  have hRe : R * ε = 2 * t - a := by rw [hR, hδ, hx]; field_simp
  have hκe : κ * ε = 1.4 * t := by rw [hκ]; field_simp
  obtain ⟨hRpos, htr⟩ := box_tail2 t ε a R κ ht1 hε0 he2 ha1 hRe hκe
  have hrg := box_tail1 t ε ht1 hε0 hε1
  have he2p : 0 < ε ^ 2 := by positivity
  refine ⟨hRpos, hD, ?_⟩
  have hfin : 0.001 * ε ^ 2 + 11.7 * ε ^ 2 * exp (-2 * x ^ 2) + 0.0001 * ε ^ 2
      < 13.6 * ε ^ 2 * exp (-2 * x ^ 2) := by
    nlinarith [mul_pos he2p hE2pos]
  linarith


noncomputable def I1 (z : ℝ) : ℝ := 2 * S z / z

/-! ## Gain step: layer cake -/

lemma integral_exp_mul_0 {β : ℝ} (hβ : 0 < β) (y : ℝ) :
    ∫ s in (0 : ℝ)..y, exp (β * s) = (exp (β * y) - 1) / β := by
  have hd : ∀ s ∈ uIcc (0 : ℝ) y, HasDerivAt (fun s => exp (β * s) / β) (exp (β * s)) s := by
    intro s _
    have h := ((hasDerivAt_id s).const_mul β).exp.div_const β
    simp only [id, mul_one] at h
    convert h using 1; field_simp
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt hd
    ((Real.continuous_exp.comp (continuous_const.mul continuous_id)).intervalIntegrable 0 y)]
  simp only [mul_zero, Real.exp_zero]; ring

/-- Layer cake: E[(e^{β min(X,R)} − 1)/β] ≤ 1/(c − β) when P(X > s) ≤ e^{−cs}. -/
lemma expect_exp_min_le {Ω : Type*} [MeasurableSpace Ω] (μ : Measure Ω) [IsProbabilityMeasure μ]
    (X : Ω → ℝ) (hXm : Measurable X) (hX0 : ∀ ω, 0 ≤ X ω) (R β c : ℝ) (hR : 0 ≤ R)
    (hβ : 0 < β) (hcβ : β < c)
    (htail : ∀ s : ℝ, 0 < s → μ {ω | s < X ω} ≤ ENNReal.ofReal (exp (-c * s))) :
    ∫ ω, (exp (β * min (X ω) R) - 1) / β ∂μ ≤ 1 / (c - β) := by
  set Y : Ω → ℝ := fun ω => min (X ω) R with hY
  have hY0 : ∀ ω, 0 ≤ Y ω := fun ω => le_min (hX0 ω) hR
  have hYm : Measurable Y := hXm.min measurable_const
  have hYX : ∀ ω, Y ω ≤ X ω := fun ω => min_le_left _ _
  have hnn : ∀ ω, 0 ≤ (exp (β * Y ω) - 1) / β := by
    intro ω; apply div_nonneg _ hβ.le
    have : 1 ≤ exp (β * Y ω) := Real.one_le_exp (mul_nonneg hβ.le (hY0 ω))
    linarith
  have hmeas : Measurable (fun ω => (exp (β * Y ω) - 1) / β) := by
    exact ((Real.measurable_exp.comp (measurable_const.mul hYm)).sub measurable_const).div_const β
  rw [integral_eq_lintegral_of_nonneg_ae (ae_of_all _ hnn) hmeas.aestronglyMeasurable]
  have hcongr : (∫⁻ ω, ENNReal.ofReal ((exp (β * Y ω) - 1) / β) ∂μ)
      = ∫⁻ ω, ENNReal.ofReal (∫ s in (0 : ℝ)..Y ω, exp (β * s)) ∂μ := by
    apply lintegral_congr; intro ω; rw [integral_exp_mul_0 hβ]
  rw [hcongr]
  have hg_int : ∀ t > (0 : ℝ), IntervalIntegrable (fun s => exp (β * s)) volume 0 t := fun t _ =>
    (Real.continuous_exp.comp (continuous_const.mul continuous_id)).intervalIntegrable 0 t
  rw [lintegral_comp_eq_lintegral_meas_lt_mul μ (ae_of_all _ hY0) hYm.aemeasurable hg_int
    (ae_of_all _ (fun s => (Real.exp_pos _).le))]
  have hbound : (∫⁻ s in Ioi (0 : ℝ), μ {a | s < Y a} * ENNReal.ofReal (exp (β * s)))
      ≤ ∫⁻ s in Ioi (0 : ℝ), ENNReal.ofReal (exp ((β - c) * s)) := by
    apply setLIntegral_mono
    · exact ENNReal.measurable_ofReal.comp (Real.measurable_exp.comp (measurable_const.mul measurable_id))
    · intro s hs
      have h1 : μ {a | s < Y a} ≤ μ {a | s < X a} :=
        measure_mono (fun a ha => lt_of_lt_of_le ha (hYX a))
      have h2 := (h1.trans (htail s hs))
      calc μ {a | s < Y a} * ENNReal.ofReal (exp (β * s))
          ≤ ENNReal.ofReal (exp (-c * s)) * ENNReal.ofReal (exp (β * s)) := by gcongr
        _ = ENNReal.ofReal (exp ((β - c) * s)) := by
          rw [← ENNReal.ofReal_mul (Real.exp_pos _).le, ← Real.exp_add]; ring_nf
  have hval : (∫⁻ s in Ioi (0 : ℝ), ENNReal.ofReal (exp ((β - c) * s)))
      = ENNReal.ofReal (1 / (c - β)) := by
    have hneg : β - c < 0 := by linarith
    rw [← ofReal_integral_eq_lintegral_ofReal (integrableOn_exp_mul_Ioi hneg 0)
      (ae_of_all _ (fun s => (Real.exp_pos _).le)), integral_exp_mul_Ioi hneg 0]
    congr 1
    rw [mul_zero, Real.exp_zero, neg_div, ← div_neg, neg_sub]
  exact ENNReal.toReal_le_of_le_ofReal (le_of_lt (one_div_pos.mpr (by linarith)))
    (hbound.trans hval.le)

/-- Pointwise gain bound: φ(r) = F x − F(x − δr) ≤ U₀/(4x)(e^{4xδ min(r,R)} − 1) + 1{r > R}. -/
lemma phi_le (F f : ℝ → ℝ) (x δ R U₀ : ℝ) (hδ : 0 < δ) (hx1 : 1 < x) (hx18 : x ≤ 1.8)
    (hRdef : R = (x - 1) / δ) (hU : 0 ≤ U₀)
    (hF01 : ∀ z, 0 ≤ F z ∧ F z ≤ 1)
    (hFup : ∀ c d : ℝ, 1 ≤ c → c ≤ d → d ≤ 1.8 → F d - F c ≤ ∫ z in c..d, f z)
    (hfint : IntervalIntegrable f volume 1 1.8)
    (hfb : ∀ z ∈ Icc (1 : ℝ) x, f z ≤ U₀ * exp (4 * x * (x - z))) (r : ℝ) (hr : 0 ≤ r) :
    F x - F (x - δ * r) ≤ U₀ / (4 * x) * (exp (4 * x * δ * min r R) - 1)
      + (if R < r then 1 else 0) := by
  have hx0 : 0 < x := by linarith
  have hδR : δ * R = x - 1 := by rw [hRdef]; field_simp
  split_ifs with hRr
  swap
  · have h : r ≤ R := not_lt.mp hRr
    rw [min_eq_left h, add_zero]
    set c := x - δ * r with hc
    have hc1 : 1 ≤ c := by
      have : δ * r ≤ δ * R := mul_le_mul_of_nonneg_left h hδ.le
      linarith
    have hcx : c ≤ x := by have := mul_nonneg hδ.le hr; linarith
    refine (hFup c x hc1 hcx hx18).trans ?_
    have hfi : IntervalIntegrable f volume c x :=
      hfint.mono_set (by
        rw [uIcc_of_le hcx, uIcc_of_le (by norm_num : (1 : ℝ) ≤ 1.8)]
        exact Icc_subset_Icc hc1 hx18)
    have hgcont : Continuous (fun z => U₀ * exp (4 * x * (x - z))) := by fun_prop
    have hmono : (∫ z in c..x, f z) ≤ ∫ z in c..x, U₀ * exp (4 * x * (x - z)) :=
      intervalIntegral.integral_mono_on hcx hfi (hgcont.intervalIntegrable c x)
        (fun z hz => hfb z ⟨hc1.trans hz.1, hz.2⟩)
    have hd : ∀ z ∈ uIcc c x, HasDerivAt (fun z => -(U₀ / (4 * x)) * exp (4 * x * (x - z)))
        (U₀ * exp (4 * x * (x - z))) z := by
      intro z _
      have h1 := (((hasDerivAt_id z).const_sub x).const_mul (4 * x)).exp.const_mul (-(U₀ / (4 * x)))
      simp only [id] at h1
      convert h1 using 1; field_simp
    rw [intervalIntegral.integral_eq_sub_of_hasDerivAt hd (hgcont.intervalIntegrable c x)] at hmono
    have hxc : x - c = δ * r := by rw [hc]; ring
    calc (∫ z in c..x, f z) ≤ _ := hmono
      _ = U₀ / (4 * x) * (exp (4 * x * δ * r) - 1) := by
        rw [hxc]; simp only [sub_self, mul_zero, Real.exp_zero]; ring_nf
  · have h := hRr
    rw [min_eq_right h.le]
    have h1 : 0 ≤ U₀ / (4 * x) * (exp (4 * x * δ * R) - 1) := by
      apply mul_nonneg (by positivity)
      have hR0 : 0 ≤ R := by rw [hRdef]; apply div_nonneg _ hδ.le; linarith
      have : 1 ≤ exp (4 * x * δ * R) := Real.one_le_exp (by positivity)
      linarith
    linarith [(hF01 x).2, (hF01 (x - δ * r)).1]

/-- Loss: F x − F(2t) ≥ (x − 2t)·m when F x − F(2t) ≥ ∫_{2t}^x f and f ≥ m on [2t, x]. -/
lemma loss_ge (F f : ℝ → ℝ) (t x m : ℝ) (h2t : 1 ≤ 2 * t) (hle : 2 * t ≤ x) (hx : x ≤ 1.8)
    (hFlo : ∫ z in (2 * t)..x, f z ≤ F x - F (2 * t))
    (hfint : IntervalIntegrable f volume 1 1.8) (hfm : ∀ z ∈ Icc (2 * t) x, m ≤ f z) :
    (x - 2 * t) * m ≤ F x - F (2 * t) := by
  refine le_trans ?_ hFlo
  have hfi : IntervalIntegrable f volume (2 * t) x :=
    hfint.mono_set (by
      rw [uIcc_of_le hle, uIcc_of_le (by norm_num : (1 : ℝ) ≤ 1.8)]
      exact Icc_subset_Icc h2t hx)
  have := intervalIntegral.integral_mono_on hle intervalIntegrable_const hfi hfm
  rw [intervalIntegral.integral_const, smul_eq_mul] at this
  exact this

/-- Upper density bound in closed form on [1, x]. -/
lemma f_up_closed (f : ℝ → ℝ) (x : ℝ) (hx18 : x ≤ 1.8)
    (hf_up : ∀ z ∈ Icc (1 : ℝ) 1.8, f z ≤ I1 z) :
    ∀ z ∈ Icc (1 : ℝ) x, f z ≤ (4 * (4 * x - 1 / x) * (1 + 0.014) * exp (-2 * x ^ 2))
      * exp (4 * x * (x - z)) := by
  intro z hz
  have hz1 := hz.1
  have hzx := hz.2
  have hz0 : 0 < z := by linarith
  refine (hf_up z ⟨hz1, hzx.trans hx18⟩).trans ?_
  unfold I1
  have hS := S_upper hz1
  have h1 : 2 * S z / z ≤ 4 * (4 * z - 1 / z) * (1 + 0.014) * exp (-2 * z ^ 2) := by
    rw [div_le_iff₀ hz0]
    have : 4 * (4 * z - 1 / z) * (1 + 0.014) * exp (-2 * z ^ 2) * z
        = 2 * (2 * (4 * z ^ 2 - 1) * exp (-2 * z ^ 2) * (1 + 0.014)) := by field_simp; ring
    rw [this]; linarith
  have h2 : 4 * z - 1 / z ≤ 4 * x - 1 / x := by
    have : 1 / x ≤ 1 / z := one_div_le_one_div_of_le hz0 hzx
    linarith
  have h3 : exp (-2 * z ^ 2) ≤ exp (-2 * x ^ 2) * exp (4 * x * (x - z)) := by
    rw [← Real.exp_add]; apply Real.exp_le_exp.mpr; nlinarith [sq_nonneg (x - z)]
  have h4 : 0 ≤ 4 * z - 1 / z := by
    have : 1 / z ≤ 1 := by rw [div_le_one hz0]; exact hz1
    linarith
  calc 2 * S z / z ≤ 4 * (4 * z - 1 / z) * (1 + 0.014) * exp (-2 * z ^ 2) := h1
    _ ≤ 4 * (4 * x - 1 / x) * (1 + 0.014) * (exp (-2 * x ^ 2) * exp (4 * x * (x - z))) := by
        have hA : 4 * (4 * z - 1 / z) * (1 + 0.014) ≤ 4 * (4 * x - 1 / x) * (1 + 0.014) := by
          linarith
        exact mul_le_mul hA h3 (Real.exp_pos _).le (by linarith)
    _ = _ := by ring

/-- Lower bound in closed form on [2t, x]: 0.935·I₁ ≥ 0.935·4(8t − 1/(2t))e^{−2x²}. -/
lemma f_lo_closed_I1 (t x : ℝ) (ht : 0.72 ≤ t) (hle : 2 * t ≤ x) :
    ∀ z ∈ Icc (2 * t) x, 4 * (8 * t - 1 / (2 * t)) * exp (-2 * x ^ 2) * (1 - 0.065)
      ≤ (1 - 0.065) * I1 z := by
  intro z hz
  have hzx := hz.2
  have hz0 : 0 < z := by linarith [hz.1]
  have hS := S_lower (show (1 : ℝ) ≤ z by linarith [hz.1])
  have hI : 4 * (4 * z - 1 / z) * exp (-2 * z ^ 2) ≤ I1 z := by
    unfold I1; rw [le_div_iff₀ hz0]
    have : 4 * (4 * z - 1 / z) * exp (-2 * z ^ 2) * z = 2 * (2 * (4 * z ^ 2 - 1) * exp (-2 * z ^ 2)) := by
      field_simp; ring
    rw [this]; linarith
  have h2 : 4 * (8 * t - 1 / (2 * t)) * exp (-2 * x ^ 2) ≤ 4 * (4 * z - 1 / z) * exp (-2 * z ^ 2) := by
    have ht0 : 0 < 2 * t := by linarith
    have ha : 8 * t - 1 / (2 * t) ≤ 4 * z - 1 / z := by
      have : 1 / z ≤ 1 / (2 * t) := one_div_le_one_div_of_le ht0 hz.1
      linarith [hz.1]
    have hb : exp (-2 * x ^ 2) ≤ exp (-2 * z ^ 2) := Real.exp_le_exp.mpr (by nlinarith)
    have hpos : 0 ≤ 8 * t - 1 / (2 * t) := by
      have : 1 / (2 * t) ≤ 1 := by rw [div_le_one ht0]; linarith
      linarith
    exact mul_le_mul (by linarith) hb (Real.exp_pos _).le (by linarith)
  calc 4 * (8 * t - 1 / (2 * t)) * exp (-2 * x ^ 2) * (1 - 0.065)
      ≤ I1 z * (1 - 0.065) := mul_le_mul_of_nonneg_right (h2.trans hI) (by norm_num)
    _ = (1 - 0.065) * I1 z := by ring

/-! ## Gain bound (DT11 §5 steps 2-3) -/

/-- E[F(x) − F(x − δT)] ≤ δU₀/(2κ − 4xδ) + e^{−2R(R+κ)} + μ(Aᶜ). -/
lemma gain_bound {Ω : Type*} [MeasurableSpace Ω] (μ : Measure Ω) [IsProbabilityMeasure μ]
    (F f : ℝ → ℝ) (x δ κ R U₀ : ℝ) (hδ0 : 0 < δ) (hx1 : 1 < x) (hx18 : x ≤ 1.8)
    (hR : R = (x - 1) / δ) (hRpos : 0 < R) (hU0 : 0 ≤ U₀) (hDpos : 0 < 2 * κ - 4 * x * δ)
    (hF01 : ∀ z, 0 ≤ F z ∧ F z ≤ 1) (hFmono : Monotone F)
    (hFup : ∀ c d : ℝ, 1 ≤ c → c ≤ d → d ≤ 1.8 → F d - F c ≤ ∫ z in c..d, f z)
    (hfint : IntervalIntegrable f volume 1 1.8)
    (hfb : ∀ z ∈ Icc (1 : ℝ) x, f z ≤ U₀ * exp (4 * x * (x - z)))
    (T X : Ω → ℝ) (A : Set Ω) (hTm : Measurable T) (hXm : Measurable X) (hAm : MeasurableSet A)
    (hX0 : ∀ ω, 0 ≤ X ω) (hTX : ∀ ω ∈ A, T ω ≤ X ω)
    (hcross : ∀ r : ℝ, 0 < r → μ {ω | r < X ω} ≤ ENNReal.ofReal (exp (-2 * r * (r + κ)))) :
    ∫ ω, (F x - F (x - δ * T ω)) ∂μ
      ≤ δ * U₀ / (2 * κ - 4 * x * δ) + exp (-2 * R * (R + κ)) + μ.real Aᶜ := by
  have hx0 : 0 < x := by linarith
  set β := 4 * x * δ with hβ
  have hβ0 : 0 < β := by positivity
  have hFm : Measurable F := hFmono.measurable
  have hφm : Measurable (fun ω => F x - F (x - δ * T ω)) :=
    measurable_const.sub (hFm.comp (measurable_const.sub (measurable_const.mul hTm)))
  have hφi : Integrable (fun ω => F x - F (x - δ * T ω)) μ :=
    Integrable.of_mem_Icc (-1) 1 hφm.aemeasurable (ae_of_all _ (fun ω =>
      ⟨by linarith [(hF01 x).1, (hF01 (x - δ * T ω)).2],
       by linarith [(hF01 x).2, (hF01 (x - δ * T ω)).1]⟩))
  -- the three pieces of the dominating function
  set P1 : Ω → ℝ := fun ω => U₀ / (4 * x) * (exp (β * min (X ω) R) - 1) with hP1
  set P2 : Ω → ℝ := fun ω => {ω | R < X ω}.indicator (1 : Ω → ℝ) ω with hP2
  set P3 : Ω → ℝ := fun ω => Aᶜ.indicator (1 : Ω → ℝ) ω with hP3
  have hpt : ∀ ω, F x - F (x - δ * T ω) ≤ P1 ω + P2 ω + P3 ω := by
    intro ω
    have hphiX := phi_le F f x δ R U₀ hδ0 hx1 hx18 hR hU0 hF01 hFup hfint hfb (X ω) (hX0 ω)
    have hind : (if R < X ω then (1 : ℝ) else 0) = P2 ω := by
      simp only [hP2, Set.indicator_apply, Set.mem_ofPred_eq, Pi.one_apply]
    rw [hind] at hphiX
    have hfirst : 0 ≤ P1 ω := by
      simp only [hP1]
      apply mul_nonneg (by positivity)
      have : 1 ≤ exp (β * min (X ω) R) :=
        Real.one_le_exp (mul_nonneg hβ0.le (le_min (hX0 ω) hRpos.le))
      linarith
    have hi1 : 0 ≤ P2 ω := Set.indicator_nonneg (fun _ _ => zero_le_one) ω
    have hP1eq : P1 ω = U₀ / (4 * x) * (exp (4 * x * δ * min (X ω) R) - 1) := rfl
    by_cases hω : ω ∈ A
    · have hmono : F (x - δ * X ω) ≤ F (x - δ * T ω) :=
        hFmono (by have := mul_le_mul_of_nonneg_left (hTX ω hω) hδ0.le; linarith)
      have hi2 : P3 ω = 0 := Set.indicator_of_notMem (by simpa using hω) _
      rw [hi2, add_zero, hP1eq]
      linarith
    · have hi2 : P3 ω = 1 := Set.indicator_of_mem (by simpa using hω) _
      rw [hi2]
      linarith [(hF01 x).2, (hF01 (x - δ * T ω)).1]
  -- integrability
  have hYm : Measurable (fun ω => (exp (β * min (X ω) R) - 1) / β) :=
    ((Real.measurable_exp.comp (measurable_const.mul (hXm.min measurable_const))).sub
      measurable_const).div_const β
  have hYi : Integrable (fun ω => (exp (β * min (X ω) R) - 1) / β) μ := by
    apply Integrable.of_mem_Icc 0 ((exp (β * R) - 1) / β) hYm.aemeasurable
    refine ae_of_all _ (fun ω => ⟨?_, ?_⟩)
    · apply div_nonneg _ hβ0.le
      have : 1 ≤ exp (β * min (X ω) R) :=
        Real.one_le_exp (mul_nonneg hβ0.le (le_min (hX0 ω) hRpos.le))
      linarith
    · apply div_le_div_of_nonneg_right _ hβ0.le
      have : exp (β * min (X ω) R) ≤ exp (β * R) :=
        Real.exp_le_exp.mpr (mul_le_mul_of_nonneg_left (min_le_right _ _) hβ0.le)
      linarith
  have hXRm : MeasurableSet {ω | R < X ω} := measurableSet_lt measurable_const hXm
  have hP1eq' : P1 = fun ω => (δ * U₀) * ((exp (β * min (X ω) R) - 1) / β) := by
    funext ω; simp only [hP1, hβ]; field_simp
  have hP1i : Integrable P1 μ := by rw [hP1eq']; exact hYi.const_mul _
  have hP2i : Integrable P2 μ := (integrable_const (1 : ℝ)).indicator hXRm
  have hP3i : Integrable P3 μ := (integrable_const (1 : ℝ)).indicator hAm.compl
  have h12i : Integrable (fun ω => P1 ω + P2 ω) μ := hP1i.add hP2i
  have h123i : Integrable (fun ω => P1 ω + P2 ω + P3 ω) μ := h12i.add hP3i
  have hint_le : ∫ ω, (F x - F (x - δ * T ω)) ∂μ ≤ ∫ ω, (P1 ω + P2 ω + P3 ω) ∂μ :=
    integral_mono hφi h123i hpt
  have hval : ∫ ω, (P1 ω + P2 ω + P3 ω) ∂μ
      = (δ * U₀) * ∫ ω, (exp (β * min (X ω) R) - 1) / β ∂μ
        + μ.real {ω | R < X ω} + μ.real Aᶜ := by
    rw [integral_add h12i hP3i, integral_add hP1i hP2i]
    have e1 : ∫ ω, P1 ω ∂μ = (δ * U₀) * ∫ ω, (exp (β * min (X ω) R) - 1) / β ∂μ := by
      rw [hP1eq', integral_const_mul]
    have e2 : ∫ ω, P2 ω ∂μ = μ.real {ω | R < X ω} := integral_indicator_one hXRm
    have e3 : ∫ ω, P3 ω ∂μ = μ.real Aᶜ := integral_indicator_one hAm.compl
    rw [e1, e2, e3]
  -- layer cake with c = 2κ
  have htail : ∀ s : ℝ, 0 < s → μ {ω | s < X ω} ≤ ENNReal.ofReal (exp (-(2 * κ) * s)) := by
    intro s hs
    refine (hcross s hs).trans (ENNReal.ofReal_le_ofReal (Real.exp_le_exp.mpr ?_))
    nlinarith [mul_nonneg hs.le hs.le]
  have hlayer := expect_exp_min_le μ X hXm hX0 R β (2 * κ) hRpos.le hβ0 (by linarith) htail
  have hδU : 0 ≤ δ * U₀ := mul_nonneg hδ0.le hU0
  have hgain : (δ * U₀) * ∫ ω, (exp (β * min (X ω) R) - 1) / β ∂μ ≤ δ * U₀ / (2 * κ - β) := by
    calc (δ * U₀) * ∫ ω, (exp (β * min (X ω) R) - 1) / β ∂μ ≤ (δ * U₀) * (1 / (2 * κ - β)) :=
          mul_le_mul_of_nonneg_left hlayer hδU
      _ = δ * U₀ / (2 * κ - β) := by ring
  have hXR : μ.real {ω | R < X ω} ≤ exp (-2 * R * (R + κ)) := by
    rw [measureReal_def]
    exact ENNReal.toReal_le_of_le_ofReal (Real.exp_pos _).le (hcross R hRpos)
  linarith

/-! ## Regularity of I₁ (used by Theorem NE and the witness) -/

lemma term_continuous (k : ℕ) : Continuous (fun z => term z k) := by unfold term; fun_prop

lemma summable_u : Summable (fun k : ℕ => 16 * ((k : ℝ) + 1) ^ 2 * exp (-2 * ((k : ℝ) + 1) ^ 2)) := by
  have hg : Summable (fun k : ℕ => 16 * exp (-1) ^ k) :=
    (summable_geometric_of_lt_one (Real.exp_pos _).le (exp_neg_lt_one (by norm_num))).mul_left 16
  refine Summable.of_nonneg_of_le (fun k => by positivity) (fun k => ?_) hg
  have hk0 : (0 : ℝ) ≤ k := Nat.cast_nonneg k
  have h1 : (k : ℝ) + 1 ≤ exp k := by have := Real.add_one_le_exp (k : ℝ); linarith
  have h2 : ((k : ℝ) + 1) ^ 2 ≤ exp (2 * k) := by
    rw [show (2 : ℝ) * k = k + k by ring, Real.exp_add]
    nlinarith [Real.exp_pos (k : ℝ)]
  have h3 : exp (2 * k) * exp (-2 * ((k : ℝ) + 1) ^ 2) ≤ exp (-1) ^ k := by
    rw [← Real.exp_nat_mul, ← Real.exp_add]; apply Real.exp_le_exp.mpr; nlinarith
  calc 16 * ((k : ℝ) + 1) ^ 2 * exp (-2 * ((k : ℝ) + 1) ^ 2)
      ≤ 16 * exp (2 * k) * exp (-2 * ((k : ℝ) + 1) ^ 2) := by gcongr
    _ = 16 * (exp (2 * k) * exp (-2 * ((k : ℝ) + 1) ^ 2)) := by ring
    _ ≤ 16 * exp (-1) ^ k := by gcongr

lemma S_continuousOn : ContinuousOn S (Icc 1 2) := by
  have h : ContinuousOn (fun z => ∑' k, term z k) (Icc 1 2) := by
    apply continuousOn_tsum (fun i => (term_continuous i).continuousOn) summable_u
    intro n z hz
    rw [Real.norm_eq_abs, abs_of_nonneg (term_nonneg hz.1 n)]
    unfold term
    have hn : (1 : ℝ) ≤ ((n : ℝ) + 1) ^ 2 := one_le_pow₀ (by have := Nat.cast_nonneg (α := ℝ) n; linarith)
    have hz1 : (1 : ℝ) ≤ z ^ 2 := one_le_pow₀ hz.1
    have hz2 : z ^ 2 ≤ 4 := by nlinarith [hz.1, hz.2]
    have he : exp (-2 * ((n : ℝ) + 1) ^ 2 * z ^ 2) ≤ exp (-2 * ((n : ℝ) + 1) ^ 2) :=
      Real.exp_le_exp.mpr (by nlinarith)
    have hp : 4 * ((n : ℝ) + 1) ^ 2 * z ^ 2 - 1 ≤ 16 * ((n : ℝ) + 1) ^ 2 := by nlinarith
    have hp0 : 0 ≤ 4 * ((n : ℝ) + 1) ^ 2 * z ^ 2 - 1 := by nlinarith
    calc (4 * ((n : ℝ) + 1) ^ 2 * z ^ 2 - 1) * exp (-2 * ((n : ℝ) + 1) ^ 2 * z ^ 2)
        ≤ (16 * ((n : ℝ) + 1) ^ 2) * exp (-2 * ((n : ℝ) + 1) ^ 2) :=
          mul_le_mul hp he (Real.exp_pos _).le (by positivity)
      _ = _ := by ring
  unfold S
  exact continuousOn_const.mul h

lemma I1_continuousOn : ContinuousOn I1 (Icc 1 2) := by
  unfold I1
  exact (continuousOn_const.mul S_continuousOn).div continuousOn_id
    (fun z hz => by have := hz.1; exact ne_of_gt (by linarith))

lemma I1_nonneg {z : ℝ} (hz : 1 ≤ z) : 0 ≤ I1 z := by
  unfold I1
  have h := S_lower hz
  have : 0 ≤ 2 * (4 * z ^ 2 - 1) * exp (-2 * z ^ 2) := by
    have : 0 ≤ 4 * z ^ 2 - 1 := by nlinarith
    positivity
  apply div_nonneg _ (by linarith); linarith

lemma Jbar_nonneg (z : ℝ) : 0 ≤ Jbar z := by
  unfold Jbar; exact mul_nonneg (by norm_num) (tsum_nonneg (fun k => (Real.exp_pos _).le))

lemma I1_intervalIntegrable {c d : ℝ} (hc : 1 ≤ c) (hcd : c ≤ d) (hd : d ≤ 2) :
    IntervalIntegrable I1 volume c d :=
  (I1_continuousOn.mono (by rw [uIcc_of_le hcd]; exact Icc_subset_Icc hc hd)).intervalIntegrable

/-! ## Theorem NE -/

/-- **Theorem NE (near-edge band, DT11).** For t ∈ [0.72, 0.88] and ε ∈ (0, 0.15], with
a = √(1 − ε²): the Lemma-B bound g of P(d_{1−ε²} ≥ t) is strictly below 1 − F(2t) = g(1, t).
Hypotheses: F is a CDF whose increments satisfy the Lemma-1 bounds UP (on [1, 1.8]) and LO
(on [1.44, 1.8]), see `lemma1_increments`; X ≥ 0 has the bridge boundary-crossing tail with κ = 1.4t/ε; T ≤ X on the measurable set A;
μ(Aᶜ), μ(E) ≤ S(0.6t/ε) (Kuiper law, G ≤ range); g obeys Lemma B. -/
theorem near_edge_band
    {Ω : Type*} [MeasurableSpace Ω] (μ : Measure Ω) [IsProbabilityMeasure μ]
    (F : ℝ → ℝ) (t ε g : ℝ)
    (ht : 0.72 ≤ t ∧ t ≤ 0.88) (hε : 0 < ε ∧ ε ≤ 0.15)
    (hF01 : ∀ z, 0 ≤ F z ∧ F z ≤ 1) (hFmono : Monotone F)
    (hUP : ∀ c d : ℝ, 1 ≤ c → c ≤ d → d ≤ 1.8 → F d - F c ≤ ∫ z in c..d, I1 z)
    (hLO : ∀ c d : ℝ, 1.44 ≤ c → c ≤ d → d ≤ 1.8 →
      (1 - 0.065) * ∫ z in c..d, I1 z ≤ F d - F c)
    (T X : Ω → ℝ) (A E : Set Ω) (hTm : Measurable T) (hXm : Measurable X)
    (hAm : MeasurableSet A)
    (hX0 : ∀ ω, 0 ≤ X ω) (hTX : ∀ ω ∈ A, T ω ≤ X ω)
    (hcross : ∀ r : ℝ, 0 < r →
      μ {ω | r < X ω} ≤ ENNReal.ofReal (exp (-2 * r * (r + 1.4 * t / ε))))
    (hA : (μ Aᶜ).toReal ≤ S (0.6 * t / ε))
    (hE : (μ E).toReal ≤ S (0.6 * t / ε))
    (hB : g ≤ (μ E).toReal +
      ∫ ω, (1 - F (2 * t / √(1 - ε ^ 2) - ε / √(1 - ε ^ 2) * T ω)) ∂μ) :
    g < 1 - F (2 * t) := by
  obtain ⟨ht1, ht2⟩ := ht
  obtain ⟨hε0, hε1⟩ := hε
  set a := √(1 - ε ^ 2) with ha
  set x := 2 * t / a with hx
  set δ := ε / a with hδ
  set κ := 1.4 * t / ε with hκ
  set R := (x - 1) / δ with hR
  obtain ⟨hRpos, hDpos, hscal⟩ := scalar_box t ε a x δ κ R ht1 ht2 hε0 hε1 ha hx hδ hκ hR
  obtain ⟨ha2, ha0, _ha1, _haL, _he2⟩ := box_a ε a hε0 hε1 ha
  have hxa : x * a = 2 * t := by rw [hx]; field_simp
  have hx2t : 2 * t ≤ x := by nlinarith
  have hx18 : x ≤ 1.8 := by
    have : x ^ 2 * a ^ 2 = 4 * t ^ 2 := by rw [← mul_pow, hxa]; ring
    nlinarith
  have hx1 : 1 < x := by linarith
  have hx0 : 0 < x := by linarith
  have hδ0 : 0 < δ := by rw [hδ]; positivity
  have hU0 : 0 ≤ 4 * (4 * x - 1 / x) * (1 + 0.014) * exp (-2 * x ^ 2) := by
    have : 0 ≤ 4 * x - 1 / x := by
      have : 1 / x ≤ 1 := by rw [div_le_one hx0]; linarith
      linarith
    positivity
  have hfint : IntervalIntegrable I1 volume 1 1.8 :=
    I1_intervalIntegrable le_rfl (by norm_num) (by norm_num)
  have hfb := f_up_closed I1 x hx18 (fun z _ => le_rfl)
  have hG := gain_bound μ F I1 x δ κ R _ hδ0 hx1 hx18 hR hRpos hU0 hDpos hF01 hFmono hUP
    hfint hfb T X A hTm hXm hAm hX0 hTX hcross
  have hfm := f_lo_closed_I1 t x ht1 hx2t
  have hFlo : ∫ z in (2 * t)..x, (1 - 0.065) * I1 z ≤ F x - F (2 * t) := by
    rw [intervalIntegral.integral_const_mul]; exact hLO (2 * t) x (by linarith) hx2t hx18
  have hloss := loss_ge F (fun z => (1 - 0.065) * I1 z) t x _ (by linarith) hx2t hx18 hFlo
    (hfint.const_mul _) hfm
  -- split the Lemma-B integral
  have hφi : Integrable (fun ω => F x - F (x - δ * T ω)) μ :=
    Integrable.of_mem_Icc (-1) 1
      (measurable_const.sub (hFmono.measurable.comp
        (measurable_const.sub (measurable_const.mul hTm)))).aemeasurable
      (ae_of_all _ (fun ω => ⟨by linarith [(hF01 x).1, (hF01 (x - δ * T ω)).2],
        by linarith [(hF01 x).2, (hF01 (x - δ * T ω)).1]⟩))
  have hsplit : ∫ ω, (1 - F (x - δ * T ω)) ∂μ = (1 - F x) + ∫ ω, (F x - F (x - δ * T ω)) ∂μ := by
    have : (fun ω => 1 - F (x - δ * T ω)) = fun ω => (1 - F x) + (F x - F (x - δ * T ω)) := by
      funext ω; ring
    rw [this, integral_add (integrable_const _) hφi, integral_const, probReal_univ, one_smul]
  have hB' : g ≤ (μ E).toReal + ((1 - F x) + ∫ ω, (F x - F (x - δ * T ω)) ∂μ) := by
    rw [← hsplit]; exact hB
  have hAc : μ.real Aᶜ ≤ S (0.6 * t / ε) := by rw [measureReal_def]; exact hA
  linarith


/-- Theorem NE in density form: a density f of F on [1, 1.8] with I₁(1 − J̄) ≤ f ≤ I₁ gives UP
and LO (pointwise Lemma 1). -/
theorem near_edge_band_density
    {Ω : Type*} [MeasurableSpace Ω] (μ : Measure Ω) [IsProbabilityMeasure μ]
    (F f : ℝ → ℝ) (t ε g : ℝ)
    (ht : 0.72 ≤ t ∧ t ≤ 0.88) (hε : 0 < ε ∧ ε ≤ 0.15)
    (hF01 : ∀ z, 0 ≤ F z ∧ F z ≤ 1) (hFmono : Monotone F)
    (hFint : ∀ c d : ℝ, 1 ≤ c → c ≤ d → d ≤ 1.8 → F d - F c = ∫ z in c..d, f z)
    (hfint : IntervalIntegrable f volume 1 1.8)
    (hf_up : ∀ z ∈ Icc (1 : ℝ) 1.8, f z ≤ I1 z)
    (hf_lo : ∀ z ∈ Icc (1.44 : ℝ) 1.8, I1 z * (1 - Jbar z) ≤ f z)
    (T X : Ω → ℝ) (A E : Set Ω) (hTm : Measurable T) (hXm : Measurable X)
    (hAm : MeasurableSet A)
    (hX0 : ∀ ω, 0 ≤ X ω) (hTX : ∀ ω ∈ A, T ω ≤ X ω)
    (hcross : ∀ r : ℝ, 0 < r →
      μ {ω | r < X ω} ≤ ENNReal.ofReal (exp (-2 * r * (r + 1.4 * t / ε))))
    (hA : (μ Aᶜ).toReal ≤ S (0.6 * t / ε))
    (hE : (μ E).toReal ≤ S (0.6 * t / ε))
    (hB : g ≤ (μ E).toReal +
      ∫ ω, (1 - F (2 * t / √(1 - ε ^ 2) - ε / √(1 - ε ^ 2) * T ω)) ∂μ) :
    g < 1 - F (2 * t) := by
  have hfi : ∀ c d : ℝ, 1 ≤ c → c ≤ d → d ≤ 1.8 → IntervalIntegrable f volume c d :=
    fun c d hc hcd hd => hfint.mono_set (by
      rw [uIcc_of_le hcd, uIcc_of_le (by norm_num : (1 : ℝ) ≤ 1.8)]
      exact Icc_subset_Icc hc hd)
  refine near_edge_band μ F t ε g ht hε hF01 hFmono (fun c d hc hcd hd => ?_)
    (fun c d hc hcd hd => ?_) T X A E hTm hXm hAm hX0 hTX hcross hA hE hB
  · rw [hFint c d hc hcd hd]
    exact intervalIntegral.integral_mono_on hcd (hfi c d hc hcd hd)
      (I1_intervalIntegrable hc hcd (by linarith))
      (fun z hz => hf_up z ⟨hc.trans hz.1, hz.2.trans hd⟩)
  · have hc1 : (1 : ℝ) ≤ c := by linarith
    rw [hFint c d hc1 hcd hd, ← intervalIntegral.integral_const_mul]
    refine intervalIntegral.integral_mono_on hcd
      ((I1_intervalIntegrable hc1 hcd (by linarith)).const_mul _) (hfi c d hc1 hcd hd)
      (fun z hz => ?_)
    have hz1 : (1.44 : ℝ) ≤ z := hc.trans hz.1
    have h1 := hf_lo z ⟨hz1, hz.2.trans hd⟩
    have h2 := Jbar_le hz1
    have h3 := I1_nonneg (show (1 : ℝ) ≤ z by linarith)
    nlinarith [mul_le_mul_of_nonneg_left h2 h3]


/-! ## Non-vacuity witness 

A concrete model in which every hypothesis of `near_edge_band` holds, with g = 1 − F(x)
(Lemma B with equality): F(z) = ∫₁^{clamp z} I₁ (so UP holds with equality and LO since
I₁ ≥ 0), a one-point probability
space, T = X = 0, A = Ω, E = ∅, t = 0.73, ε = 0.1. -/

lemma I1_integral_le_one : ∫ u in (1 : ℝ)..1.8, I1 u ≤ 1 := by
  have hb : ∀ u ∈ Icc (1 : ℝ) 1.8, I1 u ≤ 16.224 * u * exp (-2 * u ^ 2) := by
    intro u hu
    have hu0 : 0 < u := by linarith [hu.1]
    have hS := S_upper hu.1
    unfold I1; rw [div_le_iff₀ hu0]
    have : 0 ≤ exp (-2 * u ^ 2) := (Real.exp_pos _).le
    nlinarith
  have hcont : Continuous (fun u : ℝ => 16.224 * u * exp (-2 * u ^ 2)) := by fun_prop
  have h1 := intervalIntegral.integral_mono_on (by norm_num) (I1_intervalIntegrable le_rfl (by norm_num)
    (by norm_num)) (hcont.intervalIntegrable 1 1.8) hb
  have hd : ∀ u ∈ uIcc (1 : ℝ) 1.8, HasDerivAt (fun u => -4.056 * exp (-2 * u ^ 2))
      (16.224 * u * exp (-2 * u ^ 2)) u := by
    intro u _
    have h := (((hasDerivAt_pow 2 u).const_mul (-2)).exp).const_mul (-4.056)
    convert h using 1; push_cast; ring
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt hd (hcont.intervalIntegrable 1 1.8)] at h1
  have he2 : exp (-2 * (1 : ℝ) ^ 2) ≤ 1 / 7 := by
    rw [show (-2 * (1 : ℝ) ^ 2) = -(2 : ℕ) by norm_num, Real.exp_neg, inv_eq_one_div]
    apply one_div_le_one_div_of_le (by norm_num)
    have := exp_nat_ge 2; norm_num at this ⊢; linarith
  have : 0 ≤ exp (-2 * (1.8 : ℝ) ^ 2) := (Real.exp_pos _).le
  linarith

/-- The clamped CDF of the witness. -/
noncomputable def F₀ (z : ℝ) : ℝ := ∫ u in (1 : ℝ)..(max 1 (min z 1.8)), I1 u

lemma clamp_mem (z : ℝ) : 1 ≤ max 1 (min z 1.8) ∧ max 1 (min z 1.8) ≤ 1.8 :=
  ⟨le_max_left _ _, max_le (by norm_num) (min_le_right _ _)⟩

lemma F₀_diff {c d : ℝ} (hc : 1 ≤ c) (hcd : c ≤ d) (hd : d ≤ 1.8) :
    F₀ d - F₀ c = ∫ u in c..d, I1 u := by
  have ec : max 1 (min c 1.8) = c := by
    rw [min_eq_left (hcd.trans hd), max_eq_right hc]
  have ed : max 1 (min d 1.8) = d := by rw [min_eq_left hd, max_eq_right (hc.trans hcd)]
  unfold F₀; rw [ec, ed]
  exact intervalIntegral.integral_interval_sub_left
    (I1_intervalIntegrable le_rfl (hc.trans hcd) (by linarith))
    (I1_intervalIntegrable le_rfl hc (by linarith))

lemma F₀_mono : Monotone F₀ := by
  intro z z' hzz'
  have hc := clamp_mem z
  have hc' := clamp_mem z'
  have hle : max 1 (min z 1.8) ≤ max 1 (min z' 1.8) := max_le_max le_rfl (min_le_min hzz' le_rfl)
  have : F₀ z' - F₀ z = ∫ u in (max 1 (min z 1.8))..(max 1 (min z' 1.8)), I1 u := by
    unfold F₀
    exact intervalIntegral.integral_interval_sub_left
      (I1_intervalIntegrable le_rfl hc'.1 (by linarith [hc'.2]))
      (I1_intervalIntegrable le_rfl hc.1 (by linarith [hc.2]))
  have h0 : 0 ≤ ∫ u in (max 1 (min z 1.8))..(max 1 (min z' 1.8)), I1 u :=
    intervalIntegral.integral_nonneg hle (fun u hu => I1_nonneg (hc.1.trans hu.1))
  linarith

lemma F₀_mem (z : ℝ) : 0 ≤ F₀ z ∧ F₀ z ≤ 1 := by
  have hc := clamp_mem z
  constructor
  · unfold F₀
    exact intervalIntegral.integral_nonneg hc.1 (fun u hu => I1_nonneg hu.1)
  · have h18 : F₀ 1.8 = ∫ u in (1 : ℝ)..1.8, I1 u := by
      unfold F₀; rw [min_self, max_eq_right (by norm_num)]
    have := F₀_mono (show z ≤ max z 1.8 from le_max_left _ _)
    have h2 : F₀ (max z 1.8) = F₀ 1.8 := by
      unfold F₀; rw [min_eq_right (le_max_right _ _), min_self]
    linarith [I1_integral_le_one]

/-- NONVAC: every hypothesis of `near_edge_band` is satisfied by an explicit model,
with g = 1 − F(x) (Lemma B holding with equality). -/
theorem near_edge_band_nonvac :
    1 - F₀ (2 * 0.73 / √(1 - 0.1 ^ 2)) < 1 - F₀ (2 * 0.73) := by
  have hS1 : 0 ≤ S (0.6 * 0.73 / 0.1) := by
    have h := S_lower (show (1 : ℝ) ≤ 0.6 * 0.73 / 0.1 by norm_num)
    have : 0 ≤ 2 * (4 * (0.6 * 0.73 / 0.1 : ℝ) ^ 2 - 1) * exp (-2 * (0.6 * 0.73 / 0.1) ^ 2) := by
      apply mul_nonneg _ (Real.exp_pos _).le; norm_num
    linarith
  refine near_edge_band (Measure.dirac ()) F₀ 0.73 0.1 _ ⟨by norm_num, by norm_num⟩
    ⟨by norm_num, by norm_num⟩ F₀_mem F₀_mono (fun c d hc hcd hd => (F₀_diff hc hcd hd).le)
    (fun c d hc hcd hd => by
      rw [F₀_diff (by linarith) hcd hd]
      have : 0 ≤ ∫ u in c..d, I1 u :=
        intervalIntegral.integral_nonneg hcd (fun u hu => I1_nonneg (by linarith [hu.1]))
      linarith)
    (fun _ => 0) (fun _ => 0) Set.univ ∅ measurable_const measurable_const MeasurableSet.univ
    (fun _ => le_rfl) (fun _ _ => le_rfl)
    (fun r hr => by
      have : {ω : Unit | r < (fun _ => (0 : ℝ)) ω} = ∅ := by
        ext ω; simp only [Set.mem_ofPred_eq, Set.mem_empty_iff_false, iff_false, not_lt]; exact hr.le
      rw [this, measure_empty]; exact bot_le)
    (by simp only [Set.compl_univ, measure_empty, ENNReal.toReal_zero]; exact hS1)
    (by simp only [measure_empty, ENNReal.toReal_zero]; exact hS1)
    (by simp only [measure_empty, ENNReal.toReal_zero, zero_add, mul_zero, sub_zero,
          integral_const, probReal_univ, one_smul]; exact le_rfl)

/-! ## Lemma 1 (algebraic core) and Lemma 3 (path facts) -/


/-- Weierstrass product inequality: Π p ≥ 1 − Σ (1 − p) for p ∈ [0,1]. -/
lemma one_sub_sum_le_prod {ι : Type*} [DecidableEq ι] (s : Finset ι) (p : ι → ℝ)
    (hp : ∀ i, 0 ≤ p i ∧ p i ≤ 1) : 1 - ∑ i ∈ s, (1 - p i) ≤ ∏ i ∈ s, p i := by
  induction s using Finset.induction_on with
  | empty => simp
  | insert a s ha ih =>
    rw [Finset.prod_insert ha, Finset.sum_insert ha]
    have hX : 0 ≤ ∑ i ∈ s, (1 - p i) :=
      Finset.sum_nonneg (fun i _ => by linarith [(hp i).2])
    have h1 := (hp a).1
    have h2 := (hp a).2
    have : p a * (1 - ∑ i ∈ s, (1 - p i)) ≤ p a * ∏ i ∈ s, p i := mul_le_mul_of_nonneg_left ih h1
    nlinarith

theorem density_sandwich_finite {ι : Type*} [DecidableEq ι] (s : Finset ι) (h p : ι → ℝ)
    (hh : ∀ i, 0 ≤ h i) (hp : ∀ i, 0 ≤ p i ∧ p i ≤ 1) :
    ∑ j ∈ s, h j - ∑ j ∈ s, ∑ i ∈ s.erase j, h j * (1 - p i)
        ≤ ∑ j ∈ s, h j * ∏ i ∈ s.erase j, p i ∧
      ∑ j ∈ s, h j * ∏ i ∈ s.erase j, p i ≤ ∑ j ∈ s, h j := by
  constructor
  · rw [← Finset.sum_sub_distrib]
    apply Finset.sum_le_sum
    intro j _
    rw [← Finset.mul_sum]
    have := mul_le_mul_of_nonneg_left (one_sub_sum_le_prod (s.erase j) p hp) (hh j)
    linarith [show h j * (1 - ∑ i ∈ s.erase j, (1 - p i)) = h j - h j * ∑ i ∈ s.erase j, (1 - p i) by ring]
  · apply Finset.sum_le_sum
    intro j _
    have : ∏ i ∈ s.erase j, p i ≤ 1 :=
      Finset.prod_le_one₀ (fun i _ => (hp i).1) (fun i _ => (hp i).2)
    nlinarith [hh j]

/-- Lemma 3(a): G(b) ≤ range. Only the "greatest convex minorant" property of φ is used. -/
theorem gcm_gap_le_range (b φ : ℝ → ℝ) (m M : ℝ)
    (hm : ∀ z ∈ Icc (0 : ℝ) 1, m ≤ b z) (hM : ∀ z ∈ Icc (0 : ℝ) 1, b z ≤ M)
    (hφ : ∀ ψ : ℝ → ℝ, ConvexOn ℝ (Icc (0 : ℝ) 1) ψ →
      (∀ z ∈ Icc (0 : ℝ) 1, ψ z ≤ b z) → ∀ z ∈ Icc (0 : ℝ) 1, ψ z ≤ φ z) :
    ∀ z ∈ Icc (0 : ℝ) 1, b z - φ z ≤ M - m := by
  intro z hz
  have := hφ (fun _ => m) (convexOn_const m (convex_Icc 0 1)) hm z hz
  have := hM z hz
  linarith

/-- Lemma 3(b): Λ_b(y) ≤ X when the range is at most w₀ < 2y and X dominates β(d) − κd. -/
theorem Lambda_le_X (b : ℝ → ℝ) (y w₀ X : ℝ) (hy : w₀ < 2 * y) (hX0 : 0 ≤ X)
    (hw : ∀ s ∈ Icc (0 : ℝ) 1, ∀ x ∈ Icc (0 : ℝ) 1, b x - b s ≤ w₀)
    (hX : ∀ d ∈ Icc (0 : ℝ) 1, b (1 - d) - (2 * y - w₀) * d ≤ X) :
    max 0 (sSup {v | ∃ s x : ℝ, 0 ≤ s ∧ s < x ∧ x < 1 ∧
      v = b x + (1 - x) * (b x - b s - 2 * y) / (x - s)}) ≤ X := by
  apply max_le hX0
  apply Real.sSup_le _ hX0
  rintro v ⟨s, x, hs, hsx, hx1, rfl⟩
  have hden : 0 < x - s := by linarith
  have hden1 : x - s ≤ 1 := by linarith
  have hN : b x - b s - 2 * y ≤ -(2 * y - w₀) :=
    by linarith [hw s ⟨hs, by linarith⟩ x ⟨by linarith, hx1.le⟩]
  have hκ : 0 < 2 * y - w₀ := by linarith
  have h1x : 0 ≤ 1 - x := by linarith
  -- (1-x) N / (x-s) ≤ -κ (1-x)
  have hfrac : (1 - x) * (b x - b s - 2 * y) / (x - s) ≤ -(2 * y - w₀) * (1 - x) := by
    rw [div_le_iff₀ hden]
    have hneg : (1 - x) * (b x - b s - 2 * y) ≤ -(2 * y - w₀) * (1 - x) := by nlinarith
    have hk : 0 ≤ (2 * y - w₀) * (1 - x) := mul_nonneg hκ.le h1x
    nlinarith [mul_nonneg hk (show 0 ≤ 1 - (x - s) by linarith)]
  have hd := hX (1 - x) ⟨h1x, by linarith⟩
  rw [show (1 : ℝ) - (1 - x) = x by ring] at hd
  linarith



/-! ## Lemma 1 (density sandwich for M), DT11 §2

Published inputs as hypotheses: BP11 (M = max √L_j E_j, L = GEM(1), E_j iid excursion maxima),
the finite-maximum density formula, the GEM(1) first and ordered-second moment measures (P5),
and the excursion law (P2): 1 − F_E = S with density f_E. Lean proves the change of variables,
the bound J ≤ J̄, and the two density bounds (for the version f = lim f_n). -/

/-- Change of variables v = x/√ℓ: (0,1) → (x, ∞), |dv/dℓ| = x/(2ℓ√ℓ). -/
lemma cov_Ioo {x : ℝ} (hx : 0 < x) (G : ℝ → ℝ≥0∞) :
    ∫⁻ ℓ in Ioo (0 : ℝ) 1, ENNReal.ofReal (x / (2 * ℓ * √ℓ)) * G (x / √ℓ) = ∫⁻ v in Ioi x, G v := by
  have hderiv : ∀ ℓ ∈ Ioo (0 : ℝ) 1,
      HasDerivWithinAt (fun ℓ => x / √ℓ) (-(x / (2 * ℓ * √ℓ))) (Ioo 0 1) ℓ := by
    intro ℓ hℓ
    have hl : 0 < ℓ := hℓ.1
    have hs : 0 < √ℓ := Real.sqrt_pos.mpr hl
    have h1 := (Real.hasDerivAt_sqrt hl.ne').inv hs.ne'
    have h2 := h1.const_mul x
    refine (h2.congr_deriv ?_).hasDerivWithinAt.congr (fun y _ => by simp [div_eq_mul_inv])
      (by simp [div_eq_mul_inv])
    have : √ℓ ^ 2 = ℓ := Real.sq_sqrt hl.le
    field_simp
    rw [this]
  have hinj : InjOn (fun ℓ => x / √ℓ) (Ioo (0 : ℝ) 1) := by
    intro a ha b hb hab
    simp only at hab
    have hsa : 0 < √a := Real.sqrt_pos.mpr ha.1
    have hsb : 0 < √b := Real.sqrt_pos.mpr hb.1
    have : √a = √b := by
      field_simp at hab; nlinarith [hab]
    have := congrArg (· ^ 2) this
    simp only [Real.sq_sqrt ha.1.le, Real.sq_sqrt hb.1.le] at this
    exact this
  have himage : (fun ℓ => x / √ℓ) '' Ioo (0 : ℝ) 1 = Ioi x := by
    ext v; constructor
    · rintro ⟨ℓ, hℓ, rfl⟩
      have hs : 0 < √ℓ := Real.sqrt_pos.mpr hℓ.1
      have hs1 : √ℓ < 1 := by
        rw [show (1 : ℝ) = √1 by simp]; exact Real.sqrt_lt_sqrt hℓ.1.le hℓ.2
      show x < x / √ℓ
      rw [lt_div_iff₀ hs]; nlinarith
    · intro hv
      have hv0 : 0 < v := lt_trans hx hv
      refine ⟨(x / v) ^ 2, ⟨by positivity, ?_⟩, ?_⟩
      · rw [div_pow, div_lt_one (by positivity)]; exact pow_lt_pow_left₀ hv hx.le (by norm_num)
      · show x / √((x / v) ^ 2) = v
        rw [Real.sqrt_sq (by positivity)]; field_simp
  rw [← himage, lintegral_image_eq_lintegral_abs_deriv_mul measurableSet_Ioo hderiv hinj G]
  apply setLIntegral_congr_fun measurableSet_Ioo
  intro ℓ hℓ
  have hs : 0 < √ℓ := Real.sqrt_pos.mpr hℓ.1
  have hl : 0 < ℓ := hℓ.1
  simp only [abs_neg]
  rw [abs_of_pos (by positivity)]

/-- ∫_x^∞ 16c·v·e^{−2cv²} dv = 4 e^{−2cx²} (x ≥ 0, c > 0). -/
lemma lintegral_gauss_term (x c : ℝ) (hx : 0 ≤ x) (hc : 0 < c) :
    ∫⁻ v in Ioi x, ENNReal.ofReal (16 * c * v * exp (-2 * c * v ^ 2))
      = ENNReal.ofReal (4 * exp (-2 * c * x ^ 2)) := by
  have hd : ∀ v ∈ Ici x, HasDerivAt (fun v => -4 * exp (-2 * c * v ^ 2))
      (16 * c * v * exp (-2 * c * v ^ 2)) v := by
    intro v _
    have h := ((hasDerivAt_pow 2 v).const_mul (-2 * c)).exp.const_mul (-4)
    convert h using 1; push_cast; ring
  have hlim : Filter.Tendsto (fun v => -4 * exp (-2 * c * v ^ 2)) Filter.atTop (nhds 0) := by
    have h1 : Filter.Tendsto (fun v : ℝ => -2 * c * v ^ 2) Filter.atTop Filter.atBot := by
      have : Filter.Tendsto (fun v : ℝ => 2 * c * v ^ 2) Filter.atTop Filter.atTop :=
        Filter.Tendsto.const_mul_atTop (by positivity) (Filter.tendsto_pow_atTop (by norm_num))
      refine (Filter.tendsto_neg_atTop_atBot.comp this).congr (fun v => ?_)
      simp only [Function.comp]; ring
    have := (Real.tendsto_exp_atBot.comp h1).const_mul (-4)
    simpa using this
  have hpos : ∀ v ∈ Ioi x, 0 ≤ 16 * c * v * exp (-2 * c * v ^ 2) := by
    intro v hv; have : 0 ≤ v := hx.trans hv.le; positivity
  have hint := integrableOn_Ioi_deriv_of_nonneg' hd hpos hlim
  rw [← ofReal_integral_eq_lintegral_ofReal hint (ae_restrict_of_forall_mem measurableSet_Ioi hpos),
    integral_Ioi_of_hasDerivAt_of_tendsto' hd hint hlim]
  congr 1; ring


/-- The pair correction: J(x) = ∫₀¹ S(x/√m) dm/m ≤ J̄(x) for x ≥ 1. -/
lemma J_le {x : ℝ} (hx : 1 ≤ x) :
    ∫⁻ m in Ioo (0 : ℝ) 1, ENNReal.ofReal (1 / m) * ENNReal.ofReal (S (x / √m))
      ≤ ENNReal.ofReal (Jbar x) := by
  have hx0 : 0 < x := by linarith
  -- rewrite the integrand as |dv/dm| · (2 S(v)/v) at v = x/√m
  have hcongr : ∫⁻ m in Ioo (0 : ℝ) 1, ENNReal.ofReal (1 / m) * ENNReal.ofReal (S (x / √m))
      = ∫⁻ m in Ioo (0 : ℝ) 1, ENNReal.ofReal (x / (2 * m * √m))
          * ENNReal.ofReal (2 * S (x / √m) / (x / √m)) := by
    apply setLIntegral_congr_fun measurableSet_Ioo
    intro m hm
    have hm0 : 0 < m := hm.1
    have hs : 0 < √m := Real.sqrt_pos.mpr hm0
    have hS0 : 0 ≤ S (x / √m) := by
      have h1 : 1 ≤ x / √m := by
        rw [le_div_iff₀ hs]
        have : √m < 1 := by rw [show (1 : ℝ) = √1 by simp]; exact Real.sqrt_lt_sqrt hm0.le hm.2
        nlinarith
      have h := S_lower h1
      have : 0 ≤ 2 * (4 * (x / √m) ^ 2 - 1) * exp (-2 * (x / √m) ^ 2) := by
        apply mul_nonneg _ (Real.exp_pos _).le; nlinarith
      linarith
    beta_reduce
    rw [← ENNReal.ofReal_mul (by positivity), ← ENNReal.ofReal_mul (by positivity)]
    congr 1
    have : √m ^ 2 = m := Real.sq_sqrt hm0.le
    field_simp
  rw [hcongr, cov_Ioo hx0 (fun v => ENNReal.ofReal (2 * S v / v))]
  -- pointwise series bound on (x, ∞)
  set b : ℕ → ℝ → ℝ := fun k v => 16 * ((k : ℝ) + 1) ^ 2 * v * exp (-2 * ((k : ℝ) + 1) ^ 2 * v ^ 2)
    with hb
  have hbsum : ∀ v, 1 ≤ v → Summable (fun k => b k v) := by
    intro v hv
    refine Summable.of_nonneg_of_le (fun k => ?_) (fun k => ?_) (summable_u.mul_left v)
    · simp only [hb]; have : 0 ≤ v := by linarith
      positivity
    · simp only [hb]
      have hk : (1 : ℝ) ≤ ((k : ℝ) + 1) ^ 2 := one_le_pow₀ (by have := Nat.cast_nonneg (α := ℝ) k; linarith)
      have hv2 : (1 : ℝ) ≤ v ^ 2 := one_le_pow₀ hv
      have he : exp (-2 * ((k : ℝ) + 1) ^ 2 * v ^ 2) ≤ exp (-2 * ((k : ℝ) + 1) ^ 2) :=
        Real.exp_le_exp.mpr (by nlinarith)
      have : 0 ≤ 16 * ((k : ℝ) + 1) ^ 2 * v := by have : 0 ≤ v := by linarith
                                                  positivity
      calc 16 * ((k : ℝ) + 1) ^ 2 * v * exp (-2 * ((k : ℝ) + 1) ^ 2 * v ^ 2)
          ≤ 16 * ((k : ℝ) + 1) ^ 2 * v * exp (-2 * ((k : ℝ) + 1) ^ 2) := by gcongr
        _ = v * (16 * ((k : ℝ) + 1) ^ 2 * exp (-2 * ((k : ℝ) + 1) ^ 2)) := by ring
  have hpt : ∀ v ∈ Ioi x, ENNReal.ofReal (2 * S v / v) ≤ ∑' k, ENNReal.ofReal (b k v) := by
    intro v hv
    have hv1 : 1 ≤ v := hx.trans hv.le
    have hv0 : 0 < v := by linarith
    rw [← ENNReal.ofReal_tsum_of_nonneg (fun k => by simp only [hb]; positivity) (hbsum v hv1)]
    apply ENNReal.ofReal_le_ofReal
    have h1 : 2 * S v = ∑' k, 4 * term v k := by unfold S; rw [tsum_mul_left]; ring
    have h2 : ∑' k, 4 * term v k ≤ ∑' k, b k v * v := by
      refine Summable.tsum_le_tsum (fun k => ?_) ((summable_term hv1).mul_left 4)
        ((hbsum v hv1).mul_right v)
      simp only [hb]; unfold term
      have : 0 ≤ exp (-2 * ((k : ℝ) + 1) ^ 2 * v ^ 2) := (Real.exp_pos _).le
      nlinarith
    rw [div_le_iff₀ hv0, h1, ← tsum_mul_right]
    exact h2
  calc ∫⁻ v in Ioi x, ENNReal.ofReal (2 * S v / v)
      ≤ ∫⁻ v in Ioi x, ∑' k, ENNReal.ofReal (b k v) := setLIntegral_mono' measurableSet_Ioi hpt
    _ = ∑' k, ∫⁻ v in Ioi x, ENNReal.ofReal (b k v) := by
        rw [lintegral_tsum (f := fun k v => ENNReal.ofReal (b k v))
          (fun k => (ENNReal.measurable_ofReal.comp (by simp only [hb]; fun_prop)).aemeasurable)]
    _ = ∑' k : ℕ, ENNReal.ofReal (4 * exp (-2 * ((k : ℝ) + 1) ^ 2 * x ^ 2)) := by
        congr 1; funext k
        have := lintegral_gauss_term x (((k : ℝ) + 1) ^ 2) hx0.le (by positivity)
        simp only [hb]
        convert this using 3
    _ = ENNReal.ofReal (Jbar x) := by
        unfold Jbar
        rw [← tsum_mul_left, ENNReal.ofReal_tsum_of_nonneg (fun k => by positivity)]
        refine Summable.mul_left 4 (Summable.of_nonneg_of_le (fun k => (Real.exp_pos _).le)
          (fun k => ?_) (summable_u.div_const 16))
        have hk : (1 : ℝ) ≤ ((k : ℝ) + 1) ^ 2 := one_le_pow₀ (by have := Nat.cast_nonneg (α := ℝ) k; linarith)
        have hx2 : (1 : ℝ) ≤ x ^ 2 := one_le_pow₀ hx
        have he : exp (-2 * ((k : ℝ) + 1) ^ 2 * x ^ 2) ≤ exp (-2 * ((k : ℝ) + 1) ^ 2) :=
          Real.exp_le_exp.mpr (by nlinarith)
        calc exp (-2 * ((k : ℝ) + 1) ^ 2 * x ^ 2) ≤ ((k : ℝ) + 1) ^ 2 * exp (-2 * ((k : ℝ) + 1) ^ 2) := by
              nlinarith [Real.exp_pos (-2 * ((k : ℝ) + 1) ^ 2)]
          _ = 16 * ((k : ℝ) + 1) ^ 2 * exp (-2 * ((k : ℝ) + 1) ^ 2) / 16 := by ring


/-! ## Lemma 1, part A: the excursion law (P2) and the first moment (P5a) -/

lemma S_nonneg {z : ℝ} (hz : 1 ≤ z) : 0 ≤ S z := by
  have h := S_lower hz
  have : 0 ≤ 2 * (4 * z ^ 2 - 1) * exp (-2 * z ^ 2) := by
    have : 0 ≤ 4 * z ^ 2 - 1 := by nlinarith
    positivity
  linarith

lemma S_tendsto_zero : Filter.Tendsto S Filter.atTop (nhds 0) := by
  have h1 : Filter.Tendsto (fun z : ℝ => 2 * z ^ 2) Filter.atTop Filter.atTop :=
    Filter.Tendsto.const_mul_atTop (by norm_num) (Filter.tendsto_pow_atTop (by norm_num))
  have h2 := ((Real.tendsto_pow_mul_exp_neg_atTop_nhds_zero 1).comp h1).const_mul 4.056
  rw [mul_zero] at h2
  refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds h2 ?_ ?_
  · filter_upwards [Filter.eventually_ge_atTop 1] with z hz using S_nonneg hz
  · filter_upwards [Filter.eventually_ge_atTop 1] with z hz
    have h := S_upper hz
    have he : 0 ≤ exp (-2 * z ^ 2) := (Real.exp_pos _).le
    have hx : exp (-(2 * z ^ 2)) = exp (-2 * z ^ 2) := by congr 1; ring
    simp only [Function.comp, pow_one, hx]
    nlinarith

lemma pE_nonneg (P pE : ℝ → ℝ) (hPmono : Monotone P)
    (hPd : ∀ v, 1 ≤ v → HasDerivAt P (pE v) v) {v : ℝ} (hv : 1 ≤ v) : 0 ≤ pE v :=
  (hPd v hv).nonneg_of_monotone hPmono

lemma P_tendsto_one (P : ℝ → ℝ) (hPS : ∀ v, 1 ≤ v → P v = 1 - S v) :
    Filter.Tendsto P Filter.atTop (nhds 1) := by
  have h : Filter.Tendsto (fun v => 1 - S v) Filter.atTop (nhds (1 - 0)) :=
    tendsto_const_nhds.sub S_tendsto_zero
  rw [sub_zero] at h
  refine h.congr' ?_
  filter_upwards [Filter.eventually_ge_atTop 1] with v hv using (hPS v hv).symm

/-- ∫_z^∞ f_E = S(z) for z ≥ 1 (FTC + P → 1). -/
lemma tail_pE (P pE : ℝ → ℝ) (hPmono : Monotone P) (hPS : ∀ v, 1 ≤ v → P v = 1 - S v)
    (hPd : ∀ v, 1 ≤ v → HasDerivAt P (pE v) v) {z : ℝ} (hz : 1 ≤ z) :
    IntegrableOn pE (Ioi z) ∧ ∫ v in Ioi z, pE v = S z := by
  have hd : ∀ v ∈ Ici z, HasDerivAt P (pE v) v := fun v hv => hPd v (hz.trans hv)
  have hpos : ∀ v ∈ Ioi z, 0 ≤ pE v := fun v hv => pE_nonneg P pE hPmono hPd (hz.trans hv.le)
  have hlim := P_tendsto_one P hPS
  have hint := integrableOn_Ioi_deriv_of_nonneg' hd hpos hlim
  refine ⟨hint, ?_⟩
  rw [integral_Ioi_of_hasDerivAt_of_tendsto' hd hint hlim, hPS z hz]
  ring

/-- h(z, ℓ) = f_E(z/√ℓ)/√ℓ: the density at z of √ℓ·E. -/
noncomputable def hE (pE : ℝ → ℝ) (z ℓ : ℝ) : ℝ≥0∞ := ENNReal.ofReal (pE (z / √ℓ) / √ℓ)

lemma div_sqrt_ge_one {z ℓ : ℝ} (hz : 1 ≤ z) (hℓ0 : 0 < ℓ) (hℓ1 : ℓ ≤ 1) : 1 ≤ z / √ℓ := by
  have hs : 0 < √ℓ := Real.sqrt_pos.mpr hℓ0
  have hs1 : √ℓ ≤ 1 := by have := Real.sqrt_le_sqrt hℓ1; rwa [Real.sqrt_one] at this
  rw [le_div_iff₀ hs]; nlinarith

/-- First moment (P5a side): ∫₀¹ h(z, ℓ) dℓ/ℓ = I₁(z) for z ≥ 1. -/
lemma first_moment (P pE : ℝ → ℝ) (hPmono : Monotone P) (hPS : ∀ v, 1 ≤ v → P v = 1 - S v)
    (hPd : ∀ v, 1 ≤ v → HasDerivAt P (pE v) v) {z : ℝ} (hz : 1 ≤ z) :
    ∫⁻ ℓ in Ioo (0 : ℝ) 1, ENNReal.ofReal (1 / ℓ) * hE pE z ℓ = ENNReal.ofReal (I1 z) := by
  have hz0 : 0 < z := by linarith
  obtain ⟨hint, hval⟩ := tail_pE P pE hPmono hPS hPd hz
  have hcongr : ∫⁻ ℓ in Ioo (0 : ℝ) 1, ENNReal.ofReal (1 / ℓ) * hE pE z ℓ
      = ∫⁻ ℓ in Ioo (0 : ℝ) 1, ENNReal.ofReal (z / (2 * ℓ * √ℓ))
          * ENNReal.ofReal (2 * pE (z / √ℓ) / z) := by
    apply setLIntegral_congr_fun measurableSet_Ioo
    intro ℓ hℓ
    have hl : 0 < ℓ := hℓ.1
    have hs : 0 < √ℓ := Real.sqrt_pos.mpr hl
    have hp : 0 ≤ pE (z / √ℓ) := pE_nonneg P pE hPmono hPd (div_sqrt_ge_one hz hl hℓ.2.le)
    have hz0' := hz0.ne'
    have hl' := hl.ne'
    have hs' := hs.ne'
    unfold hE
    beta_reduce
    rw [← ENNReal.ofReal_mul (by positivity), ← ENNReal.ofReal_mul (by positivity)]
    congr 1
    field_simp
  rw [hcongr, cov_Ioo hz0 (fun v => ENNReal.ofReal (2 * pE v / z))]
  have hpos : ∀ v ∈ Ioi z, 0 ≤ 2 * pE v / z := fun v hv => by
    have := pE_nonneg P pE hPmono hPd (hz.trans hv.le); positivity
  have hint2 : IntegrableOn (fun v => 2 * pE v / z) (Ioi z) := (hint.const_mul 2).div_const z
  rw [← ofReal_integral_eq_lintegral_ofReal hint2
    (ae_restrict_of_forall_mem measurableSet_Ioi hpos)]
  congr 1
  rw [integral_div, integral_const_mul, hval]
  rfl

/-! ## Lemma 1, part B: one fixed stick sequence ℓ (product rule + FTC) -/

/-- Q_n(z) = Π_{i<n} P(z/√ℓ_i): the CDF of max_{i<n} √ℓ_i E_i for fixed sticks ℓ. -/
noncomputable def Qn (P : ℝ → ℝ) (ℓ : ℕ → ℝ) (n : ℕ) (z : ℝ) : ℝ :=
  ∏ i ∈ Finset.range n, P (z / √(ℓ i))

/-- h_j(z) = f_E(z/√ℓ_j)/√ℓ_j. -/
noncomputable def hj (pE : ℝ → ℝ) (ℓ : ℕ → ℝ) (j : ℕ) (z : ℝ) : ℝ := pE (z / √(ℓ j)) / √(ℓ j)

/-- f_n(z) = Σ_{j<n} h_j(z) Π_{i<n, i≠j} P(z/√ℓ_i): the derivative of Q_n. -/
noncomputable def fn (P pE : ℝ → ℝ) (ℓ : ℕ → ℝ) (n : ℕ) (z : ℝ) : ℝ :=
  ∑ j ∈ Finset.range n, hj pE ℓ j z * ∏ i ∈ (Finset.range n).erase j, P (z / √(ℓ i))

lemma Qn_hasDerivAt (P pE : ℝ → ℝ) (hPd : ∀ v, 1 ≤ v → HasDerivAt P (pE v) v)
    (ℓ : ℕ → ℝ) (hℓ : ∀ i, 0 < ℓ i ∧ ℓ i ≤ 1) (n : ℕ) {z : ℝ} (hz : 1 ≤ z) :
    HasDerivAt (Qn P ℓ n) (fn P pE ℓ n z) z := by
  have hi : ∀ i ∈ Finset.range n, HasDerivAt (fun z => P (z / √(ℓ i)))
      (pE (z / √(ℓ i)) * (1 / √(ℓ i))) z := by
    intro i _
    have h1 : HasDerivAt (fun z : ℝ => z / √(ℓ i)) (1 / √(ℓ i)) z := (hasDerivAt_id z).div_const _
    exact (hPd _ (div_sqrt_ge_one hz (hℓ i).1 (hℓ i).2)).comp z h1
  have h := HasDerivAt.fun_finsetProd hi
  unfold Qn fn
  convert h using 1
  apply Finset.sum_congr rfl
  intro j _
  unfold hj
  rw [smul_eq_mul]; ring

lemma hj_nonneg (P pE : ℝ → ℝ) (hPmono : Monotone P) (hPd : ∀ v, 1 ≤ v → HasDerivAt P (pE v) v)
    (ℓ : ℕ → ℝ) (hℓ : ∀ i, 0 < ℓ i ∧ ℓ i ≤ 1) (j : ℕ) {z : ℝ} (hz : 1 ≤ z) :
    0 ≤ hj pE ℓ j z := by
  unfold hj
  have hs : 0 < √(ℓ j) := Real.sqrt_pos.mpr (hℓ j).1
  exact div_nonneg (pE_nonneg P pE hPmono hPd (div_sqrt_ge_one hz (hℓ j).1 (hℓ j).2)) hs.le

lemma fn_nonneg (P pE : ℝ → ℝ) (hPmono : Monotone P) (hP01 : ∀ v, 0 ≤ P v ∧ P v ≤ 1)
    (hPd : ∀ v, 1 ≤ v → HasDerivAt P (pE v) v)
    (ℓ : ℕ → ℝ) (hℓ : ∀ i, 0 < ℓ i ∧ ℓ i ≤ 1) (n : ℕ) {z : ℝ} (hz : 1 ≤ z) :
    0 ≤ fn P pE ℓ n z := by
  unfold fn
  apply Finset.sum_nonneg; intro j _
  exact mul_nonneg (hj_nonneg P pE hPmono hPd ℓ hℓ j hz)
    (Finset.prod_nonneg (fun i _ => (hP01 _).1))

lemma fn_continuousOn (P pE : ℝ → ℝ) (hPd : ∀ v, 1 ≤ v → HasDerivAt P (pE v) v)
    (hpEc : ContinuousOn pE (Ici 1))
    (ℓ : ℕ → ℝ) (hℓ : ∀ i, 0 < ℓ i ∧ ℓ i ≤ 1) (n : ℕ) :
    ContinuousOn (fn P pE ℓ n) (Ici 1) := by
  have hP : ∀ i, ContinuousOn (fun z => P (z / √(ℓ i))) (Ici 1) := by
    intro i z hz
    have h1 : ContinuousAt (fun z : ℝ => z / √(ℓ i)) z := continuousAt_id.div_const _
    have h2 : ContinuousAt P (z / √(ℓ i)) :=
      (hPd _ (div_sqrt_ge_one hz (hℓ i).1 (hℓ i).2)).continuousAt
    exact (ContinuousAt.comp (g := P) (f := fun z : ℝ => z / √(ℓ i)) (x := z) h2 h1).continuousWithinAt
  have hh : ∀ j, ContinuousOn (hj pE ℓ j) (Ici 1) := by
    intro j
    unfold hj
    refine ContinuousOn.div_const ?_ _
    exact hpEc.comp (continuousOn_id.div_const _)
      (fun z hz => div_sqrt_ge_one hz (hℓ j).1 (hℓ j).2)
  unfold fn
  exact continuousOn_finsetSum _
    (fun j _ => (hh j).mul (continuousOn_finsetProd _ (fun i _ => hP i)))

/-- FTC for the finite product: Q_n(d) − Q_n(c) = ∫_c^d f_n. -/
lemma Qn_diff (P pE : ℝ → ℝ) (hPd : ∀ v, 1 ≤ v → HasDerivAt P (pE v) v)
    (hpEc : ContinuousOn pE (Ici 1))
    (ℓ : ℕ → ℝ) (hℓ : ∀ i, 0 < ℓ i ∧ ℓ i ≤ 1) (n : ℕ) {c d : ℝ} (hc : 1 ≤ c) (hcd : c ≤ d) :
    Qn P ℓ n d - Qn P ℓ n c = ∫ z in c..d, fn P pE ℓ n z := by
  have hsub : uIcc c d ⊆ Ici 1 := by
    rw [uIcc_of_le hcd]; exact fun z hz => hc.trans hz.1
  have hint : IntervalIntegrable (fn P pE ℓ n) volume c d :=
    ((fn_continuousOn P pE hPd hpEc ℓ hℓ n).mono hsub).intervalIntegrable
  exact (intervalIntegral.integral_eq_sub_of_hasDerivAt
    (fun z hz => Qn_hasDerivAt P pE hPd ℓ hℓ n (hsub hz)) hint).symm

/-! ## Lemma 1, part C: pointwise bounds f_n ≤ H_n and H_n ≤ f_n + D_n -/

/-- k(z, m) = 1 − F_E(z/√m) (= S(z/√m) when z ≥ 1, m ∈ (0, 1]). -/
noncomputable def kE (P : ℝ → ℝ) (z m : ℝ) : ℝ≥0∞ := ENNReal.ofReal (1 - P (z / √m))

lemma ofReal_fn_le_tsum (P pE : ℝ → ℝ) (hPmono : Monotone P) (hP01 : ∀ v, 0 ≤ P v ∧ P v ≤ 1)
    (hPd : ∀ v, 1 ≤ v → HasDerivAt P (pE v) v)
    (ℓ : ℕ → ℝ) (hℓ : ∀ i, 0 < ℓ i ∧ ℓ i ≤ 1) (n : ℕ) {z : ℝ} (hz : 1 ≤ z) :
    ENNReal.ofReal (fn P pE ℓ n z) ≤ ∑' j, hE pE z (ℓ j) := by
  have hh := fun j => hj_nonneg P pE hPmono hPd ℓ hℓ j hz
  have h := (density_sandwich_finite (Finset.range n) (fun j => hj pE ℓ j z)
    (fun i => P (z / √(ℓ i))) hh (fun i => hP01 _)).2
  calc ENNReal.ofReal (fn P pE ℓ n z)
      ≤ ENNReal.ofReal (∑ j ∈ Finset.range n, hj pE ℓ j z) := ENNReal.ofReal_le_ofReal h
    _ = ∑ j ∈ Finset.range n, hE pE z (ℓ j) := by
        rw [ENNReal.ofReal_sum_of_nonneg (fun j _ => hh j)]; rfl
    _ ≤ ∑' j, hE pE z (ℓ j) := ENNReal.sum_le_tsum _

lemma finite_pairs_le_tsum (a : ℕ → ℕ → ℝ≥0∞) (n : ℕ) :
    ∑ j ∈ Finset.range n, ∑ i ∈ (Finset.range n).erase j, a j i
      ≤ ∑' j, ∑' i, (if i = j then 0 else a j i) := by
  calc ∑ j ∈ Finset.range n, ∑ i ∈ (Finset.range n).erase j, a j i
      ≤ ∑ j ∈ Finset.range n, ∑' i, (if i = j then 0 else a j i) := by
        apply Finset.sum_le_sum; intro j _
        calc ∑ i ∈ (Finset.range n).erase j, a j i
            = ∑ i ∈ (Finset.range n).erase j, (if i = j then 0 else a j i) := by
              apply Finset.sum_congr rfl; intro i hi
              exact (ite_eq_right_iff.mpr fun h => absurd h (Finset.ne_of_mem_erase hi)).symm
          _ ≤ _ := ENNReal.sum_le_tsum _
    _ ≤ _ := ENNReal.sum_le_tsum _

lemma sum_le_fn_add_pairs (P pE : ℝ → ℝ) (hPmono : Monotone P) (hP01 : ∀ v, 0 ≤ P v ∧ P v ≤ 1)
    (hPd : ∀ v, 1 ≤ v → HasDerivAt P (pE v) v)
    (ℓ : ℕ → ℝ) (hℓ : ∀ i, 0 < ℓ i ∧ ℓ i ≤ 1) (n : ℕ) {z : ℝ} (hz : 1 ≤ z) :
    ∑ j ∈ Finset.range n, hE pE z (ℓ j)
      ≤ ENNReal.ofReal (fn P pE ℓ n z)
        + ∑' j, ∑' i, (if i = j then 0 else hE pE z (ℓ j) * kE P z (ℓ i)) := by
  have hh := fun j => hj_nonneg P pE hPmono hPd ℓ hℓ j hz
  have hreal := (density_sandwich_finite (Finset.range n) (fun j => hj pE ℓ j z)
    (fun i => P (z / √(ℓ i))) hh (fun i => hP01 _)).1
  have hD0 : ∀ j i, 0 ≤ hj pE ℓ j z * (1 - P (z / √(ℓ i))) := fun j i =>
    mul_nonneg (hh j) (by linarith [(hP01 (z / √(ℓ i))).2])
  have hDs : 0 ≤ ∑ j ∈ Finset.range n, ∑ i ∈ (Finset.range n).erase j,
      hj pE ℓ j z * (1 - P (z / √(ℓ i))) :=
    Finset.sum_nonneg fun j _ => Finset.sum_nonneg fun i _ => hD0 j i
  have hfn0 := fn_nonneg P pE hPmono hP01 hPd ℓ hℓ n hz
  calc ∑ j ∈ Finset.range n, hE pE z (ℓ j)
      = ENNReal.ofReal (∑ j ∈ Finset.range n, hj pE ℓ j z) := by
        rw [ENNReal.ofReal_sum_of_nonneg (fun j _ => hh j)]; rfl
    _ ≤ ENNReal.ofReal (fn P pE ℓ n z + ∑ j ∈ Finset.range n, ∑ i ∈ (Finset.range n).erase j,
          hj pE ℓ j z * (1 - P (z / √(ℓ i)))) := by
        apply ENNReal.ofReal_le_ofReal
        unfold fn
        linarith
    _ = ENNReal.ofReal (fn P pE ℓ n z) + ∑ j ∈ Finset.range n, ∑ i ∈ (Finset.range n).erase j,
          hE pE z (ℓ j) * kE P z (ℓ i) := by
        rw [ENNReal.ofReal_add hfn0 hDs,
          ENNReal.ofReal_sum_of_nonneg (fun j _ => Finset.sum_nonneg fun i _ => hD0 j i)]
        congr 1
        apply Finset.sum_congr rfl; intro j _
        rw [ENNReal.ofReal_sum_of_nonneg (fun i _ => hD0 j i)]
        apply Finset.sum_congr rfl; intro i _
        rw [ENNReal.ofReal_mul (hh j)]; rfl
    _ ≤ _ := add_le_add le_rfl (finite_pairs_le_tsum _ n)

/-! ## Lemma 1, part D: increment form of the density sandwich (DT11 §2)

Published inputs, stated as hypotheses (paper side, no Lean proof required):
* P1 (Balabdaoui–Pitman 2011, Thm 1, + independence of L and (E_j) + {M_n ≤ z} ↓ {M ≤ z}):
  F z = lim_n E Π_{i<n} F_E(z/√L_i) for z ∈ [1, 1.8].
* P2 (Kennedy 1976): F_E = P is a CDF (monotone, values in [0, 1]) with P = 1 − S on [1, ∞)
  and a continuous density pE there.
* P5 (Pitman, GEM(1) moment measures): first moment ∫₀¹ h(ℓ) dℓ/ℓ; ordered pairs
  ∬_{ℓ+m<1} h(ℓ) k(m) dℓ dm/(ℓ m).
`0 < L j ω` for every ω is a model choice: GEM(1) sticks are positive almost surely, and a
modification on one common null set changes neither P1 nor P5.
Lean proves: (UP) F d − F c ≤ ∫_c^d I₁ on [1, 1.8] and (LO) 0.935 ∫_c^d I₁ ≤ F d − F c on
[1.44, 1.8]. -/

/-- The open simplex {ℓ > 0, m > 0, ℓ + m < 1}. -/
def simplex2 : Set (ℝ × ℝ) := {p | 0 < p.1 ∧ 0 < p.2 ∧ p.1 + p.2 < 1}

lemma simplex2_sub : simplex2 ⊆ Ioo (0 : ℝ) 1 ×ˢ Ioo (0 : ℝ) 1 := by
  intro p hp
  obtain ⟨h1, h2, h3⟩ := hp
  exact ⟨⟨h1, by linarith⟩, ⟨h2, by linarith⟩⟩

/-- Ordered-pair bound: drop the simplex constraint (the integrand is ≥ 0). -/
lemma pair_le_prod (h k : ℝ → ℝ≥0∞) (hh : Measurable h) (hk : Measurable k) :
    ∫⁻ p in simplex2, ENNReal.ofReal (1 / (p.1 * p.2)) * (h p.1 * k p.2)
      ≤ (∫⁻ ℓ in Ioo (0 : ℝ) 1, ENNReal.ofReal (1 / ℓ) * h ℓ)
        * ∫⁻ m in Ioo (0 : ℝ) 1, ENNReal.ofReal (1 / m) * k m := by
  have hr : Measurable (fun x : ℝ => ENNReal.ofReal (1 / x)) :=
    (measurable_const.div measurable_id).ennreal_ofReal
  calc ∫⁻ p in simplex2, ENNReal.ofReal (1 / (p.1 * p.2)) * (h p.1 * k p.2)
      ≤ ∫⁻ p in Ioo (0 : ℝ) 1 ×ˢ Ioo (0 : ℝ) 1,
          ENNReal.ofReal (1 / (p.1 * p.2)) * (h p.1 * k p.2) := lintegral_mono_set simplex2_sub
    _ = ∫⁻ p in Ioo (0 : ℝ) 1 ×ˢ Ioo (0 : ℝ) 1,
          (ENNReal.ofReal (1 / p.1) * h p.1) * (ENNReal.ofReal (1 / p.2) * k p.2) := by
        apply setLIntegral_congr_fun (measurableSet_Ioo.prod measurableSet_Ioo)
        intro p hp
        have h1 : 0 < p.1 := hp.1.1
        have h2 : 0 < p.2 := hp.2.1
        have e : ENNReal.ofReal (1 / (p.1 * p.2))
            = ENNReal.ofReal (1 / p.1) * ENNReal.ofReal (1 / p.2) := by
          rw [← ENNReal.ofReal_mul (by positivity)]; congr 1; field_simp
        beta_reduce
        rw [e]; ring
    _ = _ := by
        rw [Measure.volume_eq_prod ℝ ℝ, ← Measure.prod_restrict]
        exact lintegral_prod_mul (f := fun ℓ => ENNReal.ofReal (1 / ℓ) * h ℓ)
          (g := fun m => ENNReal.ofReal (1 / m) * k m) (hr.mul hh).aemeasurable
          (hr.mul hk).aemeasurable

/-- **Lemma 1 (increment form).** -/
theorem lemma1_increments
    {Ω : Type*} [MeasurableSpace Ω] (ν : Measure Ω) [IsProbabilityMeasure ν]
    (L : ℕ → Ω → ℝ) (hLm : ∀ j, Measurable (L j)) (hL : ∀ j ω, 0 < L j ω ∧ L j ω ≤ 1)
    (P pE : ℝ → ℝ) (hPmono : Monotone P) (hP01 : ∀ v, 0 ≤ P v ∧ P v ≤ 1)
    (hPS : ∀ v, 1 ≤ v → P v = 1 - S v)
    (hPd : ∀ v, 1 ≤ v → HasDerivAt P (pE v) v)
    (hpEc : ContinuousOn pE (Ici 1)) (hpEm : Measurable pE)
    (hP5a : ∀ h : ℝ → ℝ≥0∞, Measurable h →
      ∫⁻ ω, ∑' j, h (L j ω) ∂ν = ∫⁻ ℓ in Ioo (0 : ℝ) 1, ENNReal.ofReal (1 / ℓ) * h ℓ)
    (hP5b : ∀ h k : ℝ → ℝ≥0∞, Measurable h → Measurable k →
      ∫⁻ ω, ∑' j, ∑' i, (if i = j then 0 else h (L j ω) * k (L i ω)) ∂ν
        = ∫⁻ p in simplex2, ENNReal.ofReal (1 / (p.1 * p.2)) * (h p.1 * k p.2))
    (F : ℝ → ℝ)
    (hP1 : ∀ z ∈ Icc (1 : ℝ) 1.8, Filter.Tendsto
      (fun n => ∫ ω, ∏ i ∈ Finset.range n, P (z / √(L i ω)) ∂ν) Filter.atTop (nhds (F z))) :
    (∀ c d : ℝ, 1 ≤ c → c ≤ d → d ≤ 1.8 → F d - F c ≤ ∫ z in c..d, I1 z) ∧
    (∀ c d : ℝ, 1.44 ≤ c → c ≤ d → d ≤ 1.8 →
      (1 - 0.065) * ∫ z in c..d, I1 z ≤ F d - F c) := by
  have hPm : Measurable P := hPmono.measurable
  have hℓ : ∀ ω i, 0 < L i ω ∧ L i ω ≤ 1 := fun ω i => hL i ω
  -- measurability
  have hq : ∀ i, Measurable (fun p : ℝ × Ω => p.1 / √(L i p.2)) := fun i =>
    measurable_fst.div ((hLm i).comp measurable_snd).sqrt
  have hfnm : ∀ n, Measurable
      (fun p : ℝ × Ω => ENNReal.ofReal (fn P pE (fun i => L i p.2) n p.1)) := by
    intro n
    refine Measurable.ennreal_ofReal ?_
    unfold fn hj
    refine Finset.measurable_sum _ (fun j _ => ?_)
    exact ((hpEm.comp (hq j)).div ((hLm j).comp measurable_snd).sqrt).mul
      (Finset.measurable_prod _ (fun i _ => hPm.comp (hq i)))
  have hhm : ∀ j, Measurable (fun p : ℝ × Ω => hE pE p.1 (L j p.2)) := by
    intro j; unfold hE
    exact ((hpEm.comp (hq j)).div ((hLm j).comp measurable_snd).sqrt).ennreal_ofReal
  have hEm : ∀ z, Measurable (hE pE z) := by
    intro z; unfold hE
    exact ((hpEm.comp (measurable_const.div measurable_id.sqrt)).div
      measurable_id.sqrt).ennreal_ofReal
  have hkm : ∀ z, Measurable (kE P z) := by
    intro z; unfold kE
    exact (measurable_const.sub (hPm.comp (measurable_const.div measurable_id.sqrt))).ennreal_ofReal
  have hQm : ∀ n z, Measurable (fun ω => Qn P (fun i => L i ω) n z) := by
    intro n z; unfold Qn
    exact Finset.measurable_prod _ (fun i _ => hPm.comp (measurable_const.div (hLm i).sqrt))
  have hQi : ∀ n z, Integrable (fun ω => Qn P (fun i => L i ω) n z) ν := by
    intro n z
    refine Integrable.of_mem_Icc 0 1 (hQm n z).aemeasurable (ae_of_all _ (fun ω => ?_))
    unfold Qn
    exact ⟨Finset.prod_nonneg (fun i _ => (hP01 _).1),
      Finset.prod_le_one₀ (fun i _ => (hP01 _).1) (fun i _ => (hP01 _).2)⟩
  -- (1) increments of G_n = E Q_n as a double lintegral (FTC per ω, then Tonelli)
  have hGlin : ∀ n c d, 1 ≤ c → c ≤ d →
      0 ≤ (∫ ω, Qn P (fun i => L i ω) n d ∂ν) - ∫ ω, Qn P (fun i => L i ω) n c ∂ν ∧
      ENNReal.ofReal ((∫ ω, Qn P (fun i => L i ω) n d ∂ν) - ∫ ω, Qn P (fun i => L i ω) n c ∂ν)
        = ∫⁻ z in Ioc c d, ∫⁻ ω, ENNReal.ofReal (fn P pE (fun i => L i ω) n z) ∂ν := by
    intro n c d hc hcd
    have hdiff : ∀ ω, Qn P (fun i => L i ω) n d - Qn P (fun i => L i ω) n c
        = ∫ z in c..d, fn P pE (fun i => L i ω) n z := fun ω =>
      Qn_diff P pE hPd hpEc (fun i => L i ω) (hℓ ω) n hc hcd
    have hnn : ∀ ω, 0 ≤ ∫ z in c..d, fn P pE (fun i => L i ω) n z := fun ω =>
      intervalIntegral.integral_nonneg hcd
        (fun z hz => fn_nonneg P pE hPmono hP01 hPd _ (hℓ ω) n (hc.trans hz.1))
    have heq : (∫ ω, Qn P (fun i => L i ω) n d ∂ν) - ∫ ω, Qn P (fun i => L i ω) n c ∂ν
        = ∫ ω, (∫ z in c..d, fn P pE (fun i => L i ω) n z) ∂ν := by
      rw [← integral_sub (hQi n d) (hQi n c)]
      exact integral_congr_ae (ae_of_all _ hdiff)
    have hint : Integrable (fun ω => ∫ z in c..d, fn P pE (fun i => L i ω) n z) ν :=
      ((hQi n d).sub (hQi n c)).congr (ae_of_all _ hdiff)
    refine ⟨by rw [heq]; exact integral_nonneg hnn, ?_⟩
    rw [heq, ofReal_integral_eq_lintegral_ofReal hint (ae_of_all _ hnn)]
    have hin : ∀ ω, ENNReal.ofReal (∫ z in c..d, fn P pE (fun i => L i ω) n z)
        = ∫⁻ z in Ioc c d, ENNReal.ofReal (fn P pE (fun i => L i ω) n z) := by
      intro ω
      have hsub : uIcc c d ⊆ Ici 1 := by
        rw [uIcc_of_le hcd]; exact fun z hz => hc.trans hz.1
      have hii : IntervalIntegrable (fn P pE (fun i => L i ω) n) volume c d :=
        ((fn_continuousOn P pE hPd hpEc _ (hℓ ω) n).mono hsub).intervalIntegrable
      rw [intervalIntegral.integral_of_le hcd, ofReal_integral_eq_lintegral_ofReal hii.1
        (ae_restrict_of_forall_mem measurableSet_Ioc
          (fun z hz => fn_nonneg P pE hPmono hP01 hPd _ (hℓ ω) n (hc.trans hz.1.le)))]
    simp_rw [hin]
    exact lintegral_lintegral_swap ((hfnm n).comp measurable_swap).aemeasurable
  -- (2) inner bounds (first moment, ordered pairs)
  have hfirst : ∀ z, 1 ≤ z → ∫⁻ ω, ∑' j, hE pE z (L j ω) ∂ν = ENNReal.ofReal (I1 z) := by
    intro z hz
    rw [hP5a _ (hEm z)]; exact first_moment P pE hPmono hPS hPd hz
  have hup_in : ∀ n z, 1 ≤ z →
      ∫⁻ ω, ENNReal.ofReal (fn P pE (fun i => L i ω) n z) ∂ν ≤ ENNReal.ofReal (I1 z) := by
    intro n z hz
    rw [← hfirst z hz]
    exact lintegral_mono (fun ω => ofReal_fn_le_tsum P pE hPmono hP01 hPd _ (hℓ ω) n hz)
  have hpair_in : ∀ z, 1.44 ≤ z →
      ∫⁻ ω, ∑' j, ∑' i, (if i = j then 0 else hE pE z (L j ω) * kE P z (L i ω)) ∂ν
        ≤ ENNReal.ofReal (I1 z) * ENNReal.ofReal 0.065 := by
    intro z hz
    have hz1 : (1 : ℝ) ≤ z := by linarith
    rw [hP5b _ _ (hEm z) (hkm z)]
    refine (pair_le_prod _ _ (hEm z) (hkm z)).trans ?_
    rw [first_moment P pE hPmono hPS hPd hz1]
    refine mul_le_mul' le_rfl ?_
    calc ∫⁻ m in Ioo (0 : ℝ) 1, ENNReal.ofReal (1 / m) * kE P z m
        = ∫⁻ m in Ioo (0 : ℝ) 1, ENNReal.ofReal (1 / m) * ENNReal.ofReal (S (z / √m)) := by
          apply setLIntegral_congr_fun measurableSet_Ioo
          intro m hm
          unfold kE
          beta_reduce
          rw [hPS _ (div_sqrt_ge_one hz1 hm.1 hm.2.le), sub_sub_cancel]
      _ ≤ ENNReal.ofReal (Jbar z) := J_le hz1
      _ ≤ ENNReal.ofReal 0.065 := ENNReal.ofReal_le_ofReal (Jbar_le hz)
  have hlow_in : ∀ n z, 1.44 ≤ z →
      ∑ j ∈ Finset.range n, ∫⁻ ω, hE pE z (L j ω) ∂ν
        ≤ ∫⁻ ω, ENNReal.ofReal (fn P pE (fun i => L i ω) n z) ∂ν
          + ENNReal.ofReal (I1 z) * ENNReal.ofReal 0.065 := by
    intro n z hz
    have hz1 : (1 : ℝ) ≤ z := by linarith
    have hfm : Measurable (fun ω => ENNReal.ofReal (fn P pE (fun i => L i ω) n z)) := by
      refine Measurable.ennreal_ofReal ?_
      unfold fn hj
      refine Finset.measurable_sum _ (fun j _ => ?_)
      exact ((hpEm.comp (measurable_const.div (hLm j).sqrt)).div (hLm j).sqrt).mul
        (Finset.measurable_prod _ (fun i _ => hPm.comp (measurable_const.div (hLm i).sqrt)))
    rw [← lintegral_finsetSum' (f := fun j ω => hE pE z (L j ω)) _
      (fun j _ => ((hEm z).comp (hLm j)).aemeasurable)]
    calc ∫⁻ ω, ∑ j ∈ Finset.range n, hE pE z (L j ω) ∂ν
        ≤ ∫⁻ ω, (ENNReal.ofReal (fn P pE (fun i => L i ω) n z)
            + ∑' j, ∑' i, (if i = j then 0 else hE pE z (L j ω) * kE P z (L i ω))) ∂ν :=
          lintegral_mono (fun ω => sum_le_fn_add_pairs P pE hPmono hP01 hPd _ (hℓ ω) n hz1)
      _ = ∫⁻ ω, ENNReal.ofReal (fn P pE (fun i => L i ω) n z) ∂ν
            + ∫⁻ ω, ∑' j, ∑' i, (if i = j then 0 else hE pE z (L j ω) * kE P z (L i ω)) ∂ν :=
          lintegral_add_left hfm _
      _ ≤ _ := add_le_add le_rfl (hpair_in z hz)
  -- (3) integrate over z and let n → ∞
  refine ⟨fun c d hc hcd hd => ?_, fun c d hc hcd hd => ?_⟩
  · have hI : IntervalIntegrable I1 volume c d := I1_intervalIntegrable hc hcd (by linarith)
    have hIint : ∫⁻ z in Ioc c d, ENNReal.ofReal (I1 z) = ENNReal.ofReal (∫ z in c..d, I1 z) := by
      rw [intervalIntegral.integral_of_le hcd, ofReal_integral_eq_lintegral_ofReal hI.1
        (ae_restrict_of_forall_mem measurableSet_Ioc (fun z hz => I1_nonneg (hc.trans hz.1.le)))]
    have hIpos : 0 ≤ ∫ z in c..d, I1 z :=
      intervalIntegral.integral_nonneg hcd (fun z hz => I1_nonneg (hc.trans hz.1))
    have hn : ∀ n, (∫ ω, Qn P (fun i => L i ω) n d ∂ν) - (∫ ω, Qn P (fun i => L i ω) n c ∂ν)
        ≤ ∫ z in c..d, I1 z := by
      intro n
      rw [← ENNReal.ofReal_le_ofReal_iff hIpos, (hGlin n c d hc hcd).2, ← hIint]
      exact setLIntegral_mono' measurableSet_Ioc (fun z hz => hup_in n z (hc.trans hz.1.le))
    have hlim := (hP1 d ⟨by linarith, hd⟩).sub (hP1 c ⟨hc, by linarith⟩)
    exact le_of_tendsto' hlim hn
  · have hc1 : (1 : ℝ) ≤ c := by linarith
    have hI : IntervalIntegrable I1 volume c d := I1_intervalIntegrable hc1 hcd (by linarith)
    have hIint : ∫⁻ z in Ioc c d, ENNReal.ofReal (I1 z) = ENNReal.ofReal (∫ z in c..d, I1 z) := by
      rw [intervalIntegral.integral_of_le hcd, ofReal_integral_eq_lintegral_ofReal hI.1
        (ae_restrict_of_forall_mem measurableSet_Ioc (fun z hz => I1_nonneg (hc1.trans hz.1.le)))]
    have hIpos : 0 ≤ ∫ z in c..d, I1 z :=
      intervalIntegral.integral_nonneg hcd (fun z hz => I1_nonneg (hc1.trans hz.1))
    have hKint : ∫⁻ z in Ioc c d, ENNReal.ofReal (I1 z) * ENNReal.ofReal 0.065
        = ENNReal.ofReal ((∫ z in c..d, I1 z) * 0.065) := by
      rw [lintegral_mul_const' _ _ ENNReal.ofReal_ne_top, hIint, ENNReal.ofReal_mul hIpos]
    have hbm : ∀ j, Measurable (fun z => ∫⁻ ω, hE pE z (L j ω) ∂ν) := fun j =>
      (hhm j).lintegral_prod_right'
    have hAm : ∀ n, Measurable
        (fun z => ∫⁻ ω, ENNReal.ofReal (fn P pE (fun i => L i ω) n z) ∂ν) := fun n =>
      (hfnm n).lintegral_prod_right'
    have hbn : ∀ n, ∑ j ∈ Finset.range n, ∫⁻ z in Ioc c d, ∫⁻ ω, hE pE z (L j ω) ∂ν
        ≤ ENNReal.ofReal ((∫ ω, Qn P (fun i => L i ω) n d ∂ν)
            - ∫ ω, Qn P (fun i => L i ω) n c ∂ν)
          + ENNReal.ofReal ((∫ z in c..d, I1 z) * 0.065) := by
      intro n
      calc ∑ j ∈ Finset.range n, ∫⁻ z in Ioc c d, ∫⁻ ω, hE pE z (L j ω) ∂ν
          = ∫⁻ z in Ioc c d, ∑ j ∈ Finset.range n, ∫⁻ ω, hE pE z (L j ω) ∂ν :=
            (lintegral_finsetSum' _ (fun j _ => (hbm j).aemeasurable)).symm
        _ ≤ ∫⁻ z in Ioc c d, (∫⁻ ω, ENNReal.ofReal (fn P pE (fun i => L i ω) n z) ∂ν
              + ENNReal.ofReal (I1 z) * ENNReal.ofReal 0.065) :=
            setLIntegral_mono' measurableSet_Ioc (fun z hz => hlow_in n z (hc.trans hz.1.le))
        _ = (∫⁻ z in Ioc c d, ∫⁻ ω, ENNReal.ofReal (fn P pE (fun i => L i ω) n z) ∂ν)
              + ∫⁻ z in Ioc c d, ENNReal.ofReal (I1 z) * ENNReal.ofReal 0.065 :=
            lintegral_add_left (hAm n) _
        _ = _ := by rw [← (hGlin n c d hc1 hcd).2, hKint]
    have hbsum : ∑' j, ∫⁻ z in Ioc c d, ∫⁻ ω, hE pE z (L j ω) ∂ν
        = ENNReal.ofReal (∫ z in c..d, I1 z) := by
      rw [← lintegral_tsum (f := fun j z => ∫⁻ ω, hE pE z (L j ω) ∂ν)
        (fun j => (hbm j).aemeasurable), ← hIint]
      apply setLIntegral_congr_fun measurableSet_Ioc
      intro z hz
      beta_reduce
      rw [← lintegral_tsum (f := fun j ω => hE pE z (L j ω))
        (fun j => ((hEm z).comp (hLm j)).aemeasurable)]
      exact hfirst z (hc1.trans hz.1.le)
    have hlim := (hP1 d ⟨by linarith, hd⟩).sub (hP1 c ⟨hc1, by linarith⟩)
    have hFnn : 0 ≤ F d - F c := ge_of_tendsto' hlim (fun n => (hGlin n c d hc1 hcd).1)
    have hR := (ENNReal.tendsto_ofReal hlim).add
      (tendsto_const_nhds (x := ENNReal.ofReal ((∫ z in c..d, I1 z) * 0.065)))
    have hkey := le_of_tendsto_of_tendsto'
      (ENNReal.tendsto_nat_tsum (fun j => ∫⁻ z in Ioc c d, ∫⁻ ω, hE pE z (L j ω) ∂ν)) hR hbn
    have hK0 : 0 ≤ (∫ z in c..d, I1 z) * 0.065 := mul_nonneg hIpos (by norm_num)
    rw [hbsum, ← ENNReal.ofReal_add hFnn hK0,
      ENNReal.ofReal_le_ofReal_iff (add_nonneg hFnn hK0)] at hkey
    linarith

/-! ## Theorem NE from the published inputs (Lemma 1 + Theorem NE) -/

/-- **Theorem NE, end to end.** Under the published inputs P1, P2, P5 for F = F_M (via
`lemma1_increments`) and the bridge inputs P3, P4 with Lemma B (as in `near_edge_band`):
for t ∈ [0.72, 0.88] and ε ∈ (0, 0.15], g < 1 − F(2t). -/
theorem near_edge_band_published
    {Ω' : Type*} [MeasurableSpace Ω'] (ν : Measure Ω') [IsProbabilityMeasure ν]
    (L : ℕ → Ω' → ℝ) (hLm : ∀ j, Measurable (L j)) (hL : ∀ j ω, 0 < L j ω ∧ L j ω ≤ 1)
    (P pE : ℝ → ℝ) (hPmono : Monotone P) (hP01 : ∀ v, 0 ≤ P v ∧ P v ≤ 1)
    (hPS : ∀ v, 1 ≤ v → P v = 1 - S v)
    (hPd : ∀ v, 1 ≤ v → HasDerivAt P (pE v) v)
    (hpEc : ContinuousOn pE (Ici 1)) (hpEm : Measurable pE)
    (hP5a : ∀ h : ℝ → ℝ≥0∞, Measurable h →
      ∫⁻ ω, ∑' j, h (L j ω) ∂ν = ∫⁻ ℓ in Ioo (0 : ℝ) 1, ENNReal.ofReal (1 / ℓ) * h ℓ)
    (hP5b : ∀ h k : ℝ → ℝ≥0∞, Measurable h → Measurable k →
      ∫⁻ ω, ∑' j, ∑' i, (if i = j then 0 else h (L j ω) * k (L i ω)) ∂ν
        = ∫⁻ p in simplex2, ENNReal.ofReal (1 / (p.1 * p.2)) * (h p.1 * k p.2))
    (F : ℝ → ℝ) (hF01 : ∀ z, 0 ≤ F z ∧ F z ≤ 1) (hFmono : Monotone F)
    (hP1 : ∀ z ∈ Icc (1 : ℝ) 1.8, Filter.Tendsto
      (fun n => ∫ ω, ∏ i ∈ Finset.range n, P (z / √(L i ω)) ∂ν) Filter.atTop (nhds (F z)))
    {Ω : Type*} [MeasurableSpace Ω] (μ : Measure Ω) [IsProbabilityMeasure μ]
    (t ε g : ℝ) (ht : 0.72 ≤ t ∧ t ≤ 0.88) (hε : 0 < ε ∧ ε ≤ 0.15)
    (T X : Ω → ℝ) (A E : Set Ω) (hTm : Measurable T) (hXm : Measurable X)
    (hAm : MeasurableSet A)
    (hX0 : ∀ ω, 0 ≤ X ω) (hTX : ∀ ω ∈ A, T ω ≤ X ω)
    (hcross : ∀ r : ℝ, 0 < r →
      μ {ω | r < X ω} ≤ ENNReal.ofReal (exp (-2 * r * (r + 1.4 * t / ε))))
    (hA : (μ Aᶜ).toReal ≤ S (0.6 * t / ε))
    (hE : (μ E).toReal ≤ S (0.6 * t / ε))
    (hB : g ≤ (μ E).toReal +
      ∫ ω, (1 - F (2 * t / √(1 - ε ^ 2) - ε / √(1 - ε ^ 2) * T ω)) ∂μ) :
    g < 1 - F (2 * t) := by
  obtain ⟨hUP, hLO⟩ := lemma1_increments ν L hLm hL P pE hPmono hP01 hPS hPd hpEc hpEm
    hP5a hP5b F hP1
  exact near_edge_band μ F t ε g ht hε hF01 hFmono hUP hLO T X A E hTm hXm hAm hX0 hTX hcross
    hA hE hB

end NearEdge
