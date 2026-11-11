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
import Mathlib.Topology.Order.Compact
import Mathlib.Topology.Piecewise
import Mathlib.Topology.Order.DenselyOrdered
import Mathlib.Topology.Order.OrderClosed
import Mathlib.Analysis.Normed.Group.Bounded
import Mathlib.MeasureTheory.Constructions.BorelSpace.Order
import Mathlib.MeasureTheory.Constructions.BorelSpace.Real
import Mathlib.Probability.Independence.Basic
import Mathlib.MeasureTheory.Integral.Bochner.Basic
/-
DT12 Step 13 (work/DT12-LINEUP.md): Lemma F(d) and Lemma B of DT10 in Lean, joined to the
near-edge band of DT11. This file = `lean/NearEdge.lean` (namespace `NearEdge`, unchanged) +
the Seam modules A1, A2, A3, B, C1, C2, C3, D (namespace `Seam`), assembled by
`work/scratch/build_seam.py`.

Main results, for t ∈ [0.72, 0.88] and ε ∈ (0, 0.15], u = 1 − ε²:
* `Seam.near_edge_seam_published`: P(t ≤ d_u(B)) < P(t ≤ d_1(B)), where
  d_u(B) = `dK u B` = inf{e ≥ 0 : some φ continuous on [0,1], convex on [0,u], concave on [u,1]
  has |B − φ| ≤ e on [0,1]} (the sup-norm distance to K_u of DT7-UNIF §0 / DT10 §0).
  Hypotheses: P1, P2, P5 (Lemma 1 inputs, as in `NearEdge.lemma1_increments`) and the bridge
  inputs (S) `hind`, (M1) `hlaw1`, (M2) `hlawB`, (K) `hKuiper`, (P4) `hcross`.
* `Seam.near_edge_seam`: the same with Lemma 1's increment bounds UP/LO as hypotheses.

What Lean checks here that was only written before (DT10 §2, §6):
* Lemma F(a)-(c), Λ0: fits, chord gap `Gch` = envelope gap (`Gch_isLUB_envGap`), endpoint
  formula `Lhat` and the endpoint spectrum (`endpoint_spectrum`), scaling, perturbation bounds.
* Lemma S(ii) + Lemma F(d): gluing (`seam_restrict`, `seam_glue`), `seamFit_iff`,
  d_u ≤ e ⇔ SeamFit(e) (`dK_le_iff`, `lemmaF_d`), strict slack (`strict_slack`),
  the edge identity d_1 = Gch/2 (`dK_one`).
* Lemma B: deterministic core (`lemmaB_det`) and probability form (`lemmaB_prob`,
  `lemmaB_real`): {t ≤ d_u} ⊂ {2t ≤ εG(b₂)} ∪ {G(b₁) ≥ (2t − εΛ(b₂))/a}, then the product law.
* Measurability on (ℝ → ℝ) with the product σ-algebra by rational-index surrogates that agree
  with the real functionals on continuous paths (`GchQ_eq`, `LamQ_eq`, `XQ_eq`, `oscQ_eq`,
  `measurableSet_dK_ge`).
* The joint (`near_edge_seam`): Lemma B → `NearEdge.near_edge_band` → `edge_ge`.

What enters as HYPOTHESES (published facts, or our written derivations; DT12 §2 table):
  (S) clipped pieces b₁ = bL(B,u), b₂ = bR(B,u) independent (Gaussian splitting, our covariance
      computation; DT10 Lemma S(i));
  (M1)/(M2) Gch(b₁), Gch(B) have CDF F (BP11 concave-majorant gap law + GCM(b) = −LCM(−b) +
      sign symmetry of the bridge);
  (K) Kuiper/Vervaat range tail of b₂; (P4) finite-horizon linear-boundary crossing of b₂;
  P1, P2, P5 as in `NearEdge`.
  Lean does not build Brownian bridge or GEM(1): identifying B, F, L with those objects (incl.
  choosing a version with every path continuous) is the paper's instantiation step.
  `#print axioms` of every declaration: [propext, Classical.choice, Quot.sound]; no `sorry`.
  Axioms are not assumptions: a false hypothesis would make the theorem vacuous, so DT12 lists
  every hypothesis with its source.
Compile: `cd ../rough-vol/lean && lake env lean <dip-test>/lean/NearEdgeSeam.lean`.
-/
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

open Set

namespace Seam

/-! ## Module A1: the convex envelope on [0, 1] (sup of affine minorants) -/

/-- Affine minorants `x ↦ p.1 * x + p.2` of `h` on `[0, 1]`. -/
def AffMin (h : ℝ → ℝ) : Set (ℝ × ℝ) := {p | ∀ y ∈ Icc (0 : ℝ) 1, p.1 * y + p.2 ≤ h y}

/-- The convex envelope of `h` on `[0, 1]`: pointwise sup of its affine minorants. -/
noncomputable def env (h : ℝ → ℝ) (x : ℝ) : ℝ :=
  sSup ((fun p : ℝ × ℝ => p.1 * x + p.2) '' AffMin h)

lemma affMin_const {h : ℝ → ℝ} {m : ℝ} (hm : ∀ y ∈ Icc (0 : ℝ) 1, m ≤ h y) :
    ((0 : ℝ), m) ∈ AffMin h := by
  intro y hy; simpa using hm y hy

lemma env_nonempty {h : ℝ → ℝ} {m : ℝ} (hm : ∀ y ∈ Icc (0 : ℝ) 1, m ≤ h y) (x : ℝ) :
    ((fun p : ℝ × ℝ => p.1 * x + p.2) '' AffMin h).Nonempty :=
  ⟨_, ⟨_, affMin_const hm, rfl⟩⟩

lemma env_bdd {h : ℝ → ℝ} {x : ℝ} (hx : x ∈ Icc (0 : ℝ) 1) :
    BddAbove ((fun p : ℝ × ℝ => p.1 * x + p.2) '' AffMin h) := by
  refine ⟨h x, ?_⟩
  rintro _ ⟨p, hp, rfl⟩
  exact hp x hx

lemma le_env {h : ℝ → ℝ} {x : ℝ} (hx : x ∈ Icc (0 : ℝ) 1) {p : ℝ × ℝ} (hp : p ∈ AffMin h) :
    p.1 * x + p.2 ≤ env h x :=
  le_csSup (env_bdd hx) ⟨p, hp, rfl⟩

lemma env_le {h : ℝ → ℝ} {m : ℝ} (hm : ∀ y ∈ Icc (0 : ℝ) 1, m ≤ h y)
    {x : ℝ} (hx : x ∈ Icc (0 : ℝ) 1) : env h x ≤ h x :=
  csSup_le (env_nonempty hm x) (by rintro _ ⟨p, hp, rfl⟩; exact hp x hx)

lemma env_ge_const {h : ℝ → ℝ} {m : ℝ} (hm : ∀ y ∈ Icc (0 : ℝ) 1, m ≤ h y)
    {x : ℝ} (hx : x ∈ Icc (0 : ℝ) 1) : m ≤ env h x := by
  have := le_env hx (affMin_const hm); simpa using this

theorem env_convex {h : ℝ → ℝ} {m : ℝ} (hm : ∀ y ∈ Icc (0 : ℝ) 1, m ≤ h y) :
    ConvexOn ℝ (Icc 0 1) (env h) := by
  refine ⟨convex_Icc 0 1, ?_⟩
  intro x hx y hy a b ha hb hab
  apply csSup_le (env_nonempty hm _)
  rintro _ ⟨p, hp, rfl⟩
  have h1 := le_env hx hp
  have h2 := le_env hy hp
  simp only [smul_eq_mul]
  have e : p.1 * (a * x + b * y) + p.2 = a * (p.1 * x + p.2) + b * (p.1 * y + p.2) := by
    linear_combination (-p.2) * hab
  rw [e]
  exact add_le_add (mul_le_mul_of_nonneg_left h1 ha) (mul_le_mul_of_nonneg_left h2 hb)

/-- Lower semicontinuity of the envelope (it is a sup of continuous functions). -/
lemma env_lsc {h : ℝ → ℝ} {m : ℝ} (hm : ∀ y ∈ Icc (0 : ℝ) 1, m ≤ h y)
    {x₀ : ℝ} (_hx₀ : x₀ ∈ Icc (0 : ℝ) 1) {η : ℝ} (hη : 0 < η) :
    ∃ ρ > 0, ∀ x ∈ Icc (0 : ℝ) 1, |x - x₀| < ρ → env h x₀ - η < env h x := by
  obtain ⟨_, ⟨p, hp, rfl⟩, hlt⟩ :=
    exists_lt_of_lt_csSup (env_nonempty hm x₀) (show env h x₀ - η < env h x₀ by linarith)
  set γ := p.1 * x₀ + p.2 - (env h x₀ - η) with hγ
  have hγ0 : 0 < γ := by simp only at hlt; linarith
  refine ⟨γ / (|p.1| + 1), by positivity, fun x hx hxx => ?_⟩
  have hk : |p.1| * |x - x₀| ≤ |p.1| * (γ / (|p.1| + 1)) :=
    mul_le_mul_of_nonneg_left hxx.le (abs_nonneg _)
  have hk2 : |p.1| * (γ / (|p.1| + 1)) < γ := by
    rw [mul_div_assoc', div_lt_iff₀ (by positivity)]; nlinarith [abs_nonneg p.1]
  have hlin : p.1 * x₀ - p.1 * x ≤ |p.1| * |x - x₀| := by
    rw [← abs_mul]
    calc p.1 * x₀ - p.1 * x = -(p.1 * (x - x₀)) := by ring
      _ ≤ |p.1 * (x - x₀)| := neg_le_abs _
  have := le_env hx hp
  linarith

lemma convex_right {φ : ℝ → ℝ} (hφ : ConvexOn ℝ (Icc 0 1) φ) {x₀ x : ℝ}
    (hx₀ : 0 ≤ x₀) (hlt : x₀ < x) (hx1 : x ≤ 1) :
    φ x ≤ φ x₀ + (x - x₀) * ((φ 1 - φ x₀) / (1 - x₀)) := by
  have hd : 0 < 1 - x₀ := by linarith
  have ha : 0 ≤ (1 - x) / (1 - x₀) := div_nonneg (by linarith) hd.le
  have hb : 0 ≤ (x - x₀) / (1 - x₀) := div_nonneg (by linarith) hd.le
  have hab : (1 - x) / (1 - x₀) + (x - x₀) / (1 - x₀) = 1 := by field_simp; ring
  have h := hφ.2 (show x₀ ∈ Icc (0 : ℝ) 1 from ⟨hx₀, by linarith⟩)
    (show (1 : ℝ) ∈ Icc (0 : ℝ) 1 from ⟨by norm_num, le_rfl⟩) ha hb hab
  simp only [smul_eq_mul] at h
  have e : (1 - x) / (1 - x₀) * x₀ + (x - x₀) / (1 - x₀) * 1 = x := by field_simp; ring
  rw [e] at h
  have e2 : (1 - x) / (1 - x₀) * φ x₀ + (x - x₀) / (1 - x₀) * φ 1
      = φ x₀ + (x - x₀) * ((φ 1 - φ x₀) / (1 - x₀)) := by field_simp; ring
  linarith

lemma convex_left {φ : ℝ → ℝ} (hφ : ConvexOn ℝ (Icc 0 1) φ) {x₀ x : ℝ}
    (hx : 0 ≤ x) (hlt : x < x₀) (hx₀ : x₀ ≤ 1) :
    φ x ≤ φ x₀ + (x₀ - x) * ((φ 0 - φ x₀) / x₀) := by
  have hd : 0 < x₀ := by linarith
  have ha : 0 ≤ (x₀ - x) / x₀ := div_nonneg (by linarith) hd.le
  have hb : 0 ≤ x / x₀ := div_nonneg hx hd.le
  have hab : (x₀ - x) / x₀ + x / x₀ = 1 := by field_simp; ring
  have h := hφ.2 (show (0 : ℝ) ∈ Icc (0 : ℝ) 1 from ⟨le_rfl, by norm_num⟩)
    (show x₀ ∈ Icc (0 : ℝ) 1 from ⟨hd.le, hx₀⟩) ha hb hab
  simp only [smul_eq_mul] at h
  have e : (x₀ - x) / x₀ * 0 + x / x₀ * x₀ = x := by field_simp; ring
  rw [e] at h
  have e2 : (x₀ - x) / x₀ * φ 0 + x / x₀ * φ x₀ = φ x₀ + (x₀ - x) * ((φ 0 - φ x₀) / x₀) := by
    field_simp; ring
  linarith

/-- Upper semicontinuity of a convex function on [0, 1] (chord bounds toward 0 and 1). -/
lemma convex_usc {φ : ℝ → ℝ} (hφ : ConvexOn ℝ (Icc 0 1) φ) {x₀ : ℝ}
    (hx₀ : x₀ ∈ Icc (0 : ℝ) 1) {η : ℝ} (hη : 0 < η) :
    ∃ ρ > 0, ∀ x ∈ Icc (0 : ℝ) 1, |x - x₀| < ρ → φ x < φ x₀ + η := by
  set A := |(φ 1 - φ x₀) / (1 - x₀)| + |(φ 0 - φ x₀) / x₀| + 1 with hA
  have hA0 : 0 < A := by positivity
  refine ⟨η / A, by positivity, fun x hx hxx => ?_⟩
  have hρA : |x - x₀| * A < η := by
    have := mul_lt_mul_of_pos_right hxx hA0
    rwa [div_mul_cancel₀ η hA0.ne'] at this
  rcases lt_trichotomy x x₀ with hlt | heq | hgt
  · have h := convex_left hφ hx.1 hlt hx₀.2
    have h1 : (x₀ - x) * ((φ 0 - φ x₀) / x₀) ≤ |x - x₀| * A := by
      rw [abs_sub_comm, abs_of_pos (by linarith : 0 < x₀ - x)]
      apply mul_le_mul_of_nonneg_left _ (by linarith)
      have := le_abs_self ((φ 0 - φ x₀) / x₀)
      have := abs_nonneg ((φ 1 - φ x₀) / (1 - x₀)); linarith
    linarith
  · subst heq; linarith
  · have h := convex_right hφ hx₀.1 hgt hx.2
    have h1 : (x - x₀) * ((φ 1 - φ x₀) / (1 - x₀)) ≤ |x - x₀| * A := by
      rw [abs_of_pos (by linarith : 0 < x - x₀)]
      apply mul_le_mul_of_nonneg_left _ (by linarith)
      have := le_abs_self ((φ 1 - φ x₀) / (1 - x₀))
      have := abs_nonneg ((φ 0 - φ x₀) / x₀); linarith
    linarith

/-- A convex function on [0,1] that is lower semicontinuous (in the ε-ρ form) is continuous. -/
lemma continuousOn_of_lsc_convex {φ : ℝ → ℝ} (hφ : ConvexOn ℝ (Icc 0 1) φ)
    (hl : ∀ x₀ ∈ Icc (0 : ℝ) 1, ∀ η > 0, ∃ ρ > 0, ∀ x ∈ Icc (0 : ℝ) 1, |x - x₀| < ρ →
      φ x₀ - η < φ x) :
    ContinuousOn φ (Icc 0 1) := by
  intro x₀ hx₀
  rw [Metric.continuousWithinAt_iff]
  intro η hη
  obtain ⟨ρ₁, hρ₁, h1⟩ := hl x₀ hx₀ η hη
  obtain ⟨ρ₂, hρ₂, h2⟩ := convex_usc hφ hx₀ hη
  refine ⟨min ρ₁ ρ₂, lt_min hρ₁ hρ₂, fun x hx hd => ?_⟩
  rw [Real.dist_eq] at hd ⊢
  have a1 := h1 x hx (lt_of_lt_of_le hd (min_le_left _ _))
  have a2 := h2 x hx (lt_of_lt_of_le hd (min_le_right _ _))
  rw [abs_sub_lt_iff]; constructor <;> linarith

theorem env_continuousOn {h : ℝ → ℝ} {m : ℝ} (hm : ∀ y ∈ Icc (0 : ℝ) 1, m ≤ h y) :
    ContinuousOn (env h) (Icc 0 1) :=
  continuousOn_of_lsc_convex (env_convex hm) (fun _ hx₀ _ hη => env_lsc hm hx₀ hη)

/-- Separation: if every chord of `h` across `x₀ ∈ (0,1)` lies above `c` at `x₀`, and
`c ≤ h x₀`, then some affine minorant passes through `(x₀, c)`; so `c ≤ env h x₀`. -/
theorem le_env_of_chords {h : ℝ → ℝ} {x₀ c : ℝ} (hx₀ : x₀ ∈ Ioo (0 : ℝ) 1) (hc : c ≤ h x₀)
    (hch : ∀ s x', 0 ≤ s → s < x₀ → x₀ < x' → x' ≤ 1 →
      c * (x' - s) ≤ (x' - x₀) * h s + (x₀ - s) * h x') :
    c ≤ env h x₀ := by
  obtain ⟨h0, h1⟩ := hx₀
  set Ls : Set ℝ := {v | ∃ s, 0 ≤ s ∧ s < x₀ ∧ v = (c - h s) / (x₀ - s)} with hLs
  have hLne : Ls.Nonempty := ⟨_, 0, le_rfl, h0, rfl⟩
  have hRle : ∀ v ∈ Ls, ∀ x', x₀ < x' → x' ≤ 1 → v ≤ (h x' - c) / (x' - x₀) := by
    rintro v ⟨s, hs0, hs, rfl⟩ x' hx' hx'1
    have hd1 : 0 < x₀ - s := by linarith
    have hd2 : 0 < x' - x₀ := by linarith
    rw [div_le_div_iff₀ hd1 hd2]
    have := hch s x' hs0 hs hx' hx'1
    nlinarith
  have hLbdd : BddAbove Ls := ⟨(h 1 - c) / (1 - x₀), fun v hv => hRle v hv 1 h1 le_rfl⟩
  set μ := sSup Ls with hμ
  have hp : (μ, c - μ * x₀) ∈ AffMin h := by
    intro y hy
    simp only
    rcases lt_trichotomy y x₀ with hlt | heq | hgt
    · have hv : (c - h y) / (x₀ - y) ≤ μ := le_csSup hLbdd ⟨y, hy.1, hlt, rfl⟩
      have hd : 0 < x₀ - y := by linarith
      rw [div_le_iff₀ hd] at hv
      nlinarith
    · subst heq; linarith
    · have hv : μ ≤ (h y - c) / (y - x₀) := csSup_le hLne (fun v hv => hRle v hv y hgt hy.2)
      have hd : 0 < y - x₀ := by linarith
      rw [le_div_iff₀ hd] at hv
      nlinarith
  have := le_env ⟨h0.le, h1.le⟩ hp
  simp only at this
  linarith

end Seam

namespace Seam

/-! ## Module A2: fits, chord gap `Gch`, endpoint formula `Lhat` (DT10 Lemma F (a)-(c), Λ0) -/

/-- `φ` is a `y`-fit of `b`: continuous and convex on `[0,1]`, within `y` of `b` there. -/
def IsFit (b : ℝ → ℝ) (y : ℝ) (φ : ℝ → ℝ) : Prop :=
  ContinuousOn φ (Icc 0 1) ∧ ConvexOn ℝ (Icc 0 1) φ ∧ ∀ x ∈ Icc (0 : ℝ) 1, |φ x - b x| ≤ y

/-- The chord of `b` between `s` and `x'`, evaluated at `x`. -/
noncomputable def chord (b : ℝ → ℝ) (s x' x : ℝ) : ℝ :=
  (x' - x) / (x' - s) * b s + (x - s) / (x' - s) * b x'

/-- All chord defects `b x − chord`, `0 ≤ s < x < x' ≤ 1`. -/
def gapSet (b : ℝ → ℝ) : Set ℝ :=
  {v | ∃ s x x' : ℝ, 0 ≤ s ∧ s < x ∧ x < x' ∧ x' ≤ 1 ∧ v = b x - chord b s x' x}

/-- Chord gap `G(b)` (= sup(b − GCM b), proved below as `Gch_eq_envGap`). -/
noncomputable def Gch (b : ℝ → ℝ) : ℝ := max 0 (sSup (gapSet b))

/-- One term of the endpoint formula of DT10 Lemma F(c). -/
noncomputable def endTerm (b : ℝ → ℝ) (y s x : ℝ) : ℝ :=
  ((1 - s) * (b x - y) - (1 - x) * (b s + y)) / (x - s)

def endSet (b : ℝ → ℝ) (y : ℝ) : Set ℝ :=
  {v | ∃ s x : ℝ, 0 ≤ s ∧ s < x ∧ x < 1 ∧ v = endTerm b y s x}

/-- Endpoint formula `L̂_b(y)` (DT10 Lemma F(c)); proved below to be the least endpoint value
of a `y`-fit (`Lhat_le_fit`, `exists_fit_endpoint_eq`). -/
noncomputable def Lhat (b : ℝ → ℝ) (y : ℝ) : ℝ := max (b 1 - y) (sSup (endSet b y))

/-- Endpoint excess `Λ_b(y) = L̂_b(y) + y` (DT10 Lemma Λ0). -/
noncomputable def Lam (b : ℝ → ℝ) (y : ℝ) : ℝ := Lhat b y + y

lemma bound_of_cont {b : ℝ → ℝ} (hb : ContinuousOn b (Icc 0 1)) :
    ∃ C, ∀ x ∈ Icc (0 : ℝ) 1, |b x| ≤ C := by
  obtain ⟨C, hC⟩ := isCompact_Icc.exists_bound_of_continuousOn hb
  exact ⟨C, fun x hx => by simpa [Real.norm_eq_abs] using hC x hx⟩

lemma chord_weights {s x x' : ℝ} (hsx : s < x) (hxx' : x < x') :
    0 ≤ (x' - x) / (x' - s) ∧ 0 ≤ (x - s) / (x' - s) ∧
      (x' - x) / (x' - s) + (x - s) / (x' - s) = 1 := by
  have hd : 0 < x' - s := by linarith
  refine ⟨div_nonneg (by linarith) hd.le, div_nonneg (by linarith) hd.le, ?_⟩
  field_simp; ring

/-- A convex function lies below its chords. -/
lemma convex_le_chord {φ : ℝ → ℝ} (hφ : ConvexOn ℝ (Icc 0 1) φ) {s x x' : ℝ}
    (hs : 0 ≤ s) (hsx : s < x) (hxx' : x < x') (hx' : x' ≤ 1) : φ x ≤ chord φ s x' x := by
  obtain ⟨ha, hb, hab⟩ := chord_weights hsx hxx'
  have h := hφ.2 (show s ∈ Icc (0 : ℝ) 1 from ⟨hs, by linarith⟩)
    (show x' ∈ Icc (0 : ℝ) 1 from ⟨by linarith, hx'⟩) ha hb hab
  simp only [smul_eq_mul] at h
  have hd : 0 < x' - s := by linarith
  have e : (x' - x) / (x' - s) * s + (x - s) / (x' - s) * x' = x := by field_simp; ring
  rw [e] at h
  exact h

lemma chord_mono {f g : ℝ → ℝ} {s x x' : ℝ} (hsx : s < x) (hxx' : x < x')
    (hs : f s ≤ g s) (hx' : f x' ≤ g x') : chord f s x' x ≤ chord g s x' x := by
  obtain ⟨ha, hb, _⟩ := chord_weights hsx hxx'
  unfold chord
  exact add_le_add (mul_le_mul_of_nonneg_left hs ha) (mul_le_mul_of_nonneg_left hx' hb)

lemma chord_add_const (f : ℝ → ℝ) (c : ℝ) {s x x' : ℝ} (hsx : s < x) (hxx' : x < x') :
    chord (fun z => f z + c) s x' x = chord f s x' x + c := by
  obtain ⟨_, _, hab⟩ := chord_weights hsx hxx'
  unfold chord; linear_combination c * hab

lemma chord_abs_le {b : ℝ → ℝ} {C : ℝ} {s x x' : ℝ} (hsx : s < x) (hxx' : x < x')
    (hs : |b s| ≤ C) (hx' : |b x'| ≤ C) : |chord b s x' x| ≤ C := by
  obtain ⟨ha, hb, hab⟩ := chord_weights hsx hxx'
  unfold chord
  calc |(x' - x) / (x' - s) * b s + (x - s) / (x' - s) * b x'|
      ≤ |(x' - x) / (x' - s) * b s| + |(x - s) / (x' - s) * b x'| := abs_add_le _ _
    _ = (x' - x) / (x' - s) * |b s| + (x - s) / (x' - s) * |b x'| := by
        rw [abs_mul, abs_mul, abs_of_nonneg ha, abs_of_nonneg hb]
    _ ≤ (x' - x) / (x' - s) * C + (x - s) / (x' - s) * C :=
        add_le_add (mul_le_mul_of_nonneg_left hs ha) (mul_le_mul_of_nonneg_left hx' hb)
    _ = C := by linear_combination C * hab

lemma gapSet_nonempty (b : ℝ → ℝ) : (gapSet b).Nonempty :=
  ⟨_, 0, 1 / 2, 1, le_rfl, by norm_num, by norm_num, le_rfl, rfl⟩

lemma gapSet_bdd {b : ℝ → ℝ} {C : ℝ} (hC : ∀ x ∈ Icc (0 : ℝ) 1, |b x| ≤ C) :
    BddAbove (gapSet b) := by
  refine ⟨2 * C, ?_⟩
  rintro v ⟨s, x, x', hs, hsx, hxx', hx', rfl⟩
  have h1 := chord_abs_le hsx hxx' (hC s ⟨hs, by linarith⟩) (hC x' ⟨by linarith, hx'⟩)
  have h2 := hC x ⟨by linarith, by linarith⟩
  have := neg_abs_le (chord b s x' x)
  have := le_abs_self (b x)
  linarith

lemma Gch_nonneg (b : ℝ → ℝ) : 0 ≤ Gch b := le_max_left _ _

lemma gap_le_Gch {b : ℝ → ℝ} (hb : ContinuousOn b (Icc 0 1)) {s x x' : ℝ}
    (hs : 0 ≤ s) (hsx : s < x) (hxx' : x < x') (hx' : x' ≤ 1) :
    b x - chord b s x' x ≤ Gch b := by
  obtain ⟨C, hC⟩ := bound_of_cont hb
  exact (le_csSup (gapSet_bdd hC) ⟨s, x, x', hs, hsx, hxx', hx', rfl⟩).trans (le_max_right _ _)

/-- Lemma F(a), necessity: a `y`-fit forces `G(b) ≤ 2y`. -/
theorem Gch_le_of_fit {b φ : ℝ → ℝ} {y : ℝ} (hf : IsFit b y φ) : Gch b ≤ 2 * y := by
  obtain ⟨_, hcv, hab⟩ := hf
  have hy : 0 ≤ y := (abs_nonneg _).trans (hab 0 ⟨le_rfl, zero_le_one⟩)
  apply max_le (by linarith)
  apply csSup_le (gapSet_nonempty b)
  rintro v ⟨s, x, x', hs, hsx, hxx', hx', rfl⟩
  have hxI : x ∈ Icc (0 : ℝ) 1 := ⟨by linarith, by linarith⟩
  have h1 := abs_le.mp (hab x hxI)
  have h2 := convex_le_chord hcv hs hsx hxx' hx'
  have hs' := abs_le.mp (hab s ⟨hs, by linarith⟩)
  have hx'' := abs_le.mp (hab x' ⟨by linarith, hx'⟩)
  have h3 : chord φ s x' x ≤ chord (fun z => b z + y) s x' x :=
    chord_mono hsx hxx' (by show φ s ≤ b s + y; linarith) (by show φ x' ≤ b x' + y; linarith)
  rw [chord_add_const b y hsx hxx'] at h3
  linarith

/-- If `f ≤ g` on `(0,1)` and both are continuous on `[0,1]`, then `f ≤ g` on `[0,1]`. -/
lemma le_on_Icc_of_Ioo {f g : ℝ → ℝ} (hf : ContinuousOn f (Icc 0 1))
    (hg : ContinuousOn g (Icc 0 1)) (h : ∀ x ∈ Ioo (0 : ℝ) 1, f x ≤ g x) :
    ∀ x ∈ Icc (0 : ℝ) 1, f x ≤ g x := by
  have hc : closure (Ioo (0 : ℝ) 1) = Icc 0 1 := closure_Ioo (by norm_num)
  intro x hx
  rw [← hc] at hf hg hx
  exact le_on_closure h hf hg hx

/-- The chord condition of `le_env_of_chords` for the obstacle `b + y`, from `G(b) ≤ 2y`. -/
lemma chord_cond_of_gap {b : ℝ → ℝ} (hb : ContinuousOn b (Icc 0 1)) {y : ℝ} (hG : Gch b ≤ 2 * y)
    {s x₀ x' : ℝ} (hs : 0 ≤ s) (hsx : s < x₀) (hxx' : x₀ < x') (hx' : x' ≤ 1) :
    (b x₀ - y) * (x' - s) ≤ (x' - x₀) * (b s + y) + (x₀ - s) * (b x' + y) := by
  have hg := (gap_le_Gch hb hs hsx hxx' hx').trans hG
  have hd : 0 < x' - s := by linarith
  unfold chord at hg
  have e : (x' - x₀) / (x' - s) * b s + (x₀ - s) / (x' - s) * b x'
      = ((x' - x₀) * b s + (x₀ - s) * b x') / (x' - s) := by field_simp
  rw [e] at hg
  have hq : b x₀ - 2 * y ≤ ((x' - x₀) * b s + (x₀ - s) * b x') / (x' - s) := by linarith
  rw [le_div_iff₀ hd] at hq
  nlinarith

/-- Lemma F(a), sufficiency: if `G(b) ≤ 2y` then `env(b + y)` is a `y`-fit (the maximal fit). -/
theorem fit_env {b : ℝ → ℝ} (hb : ContinuousOn b (Icc 0 1)) {y : ℝ} (hy : 0 ≤ y)
    (hG : Gch b ≤ 2 * y) : IsFit b y (env (fun z => b z + y)) := by
  obtain ⟨C, hC⟩ := bound_of_cont hb
  have hm : ∀ z ∈ Icc (0 : ℝ) 1, -C + y ≤ b z + y := fun z hz => by
    have := neg_abs_le (b z); have := hC z hz; linarith
  have hcont := env_continuousOn hm
  refine ⟨hcont, env_convex hm, ?_⟩
  have hlow : ∀ x ∈ Icc (0 : ℝ) 1, b x - y ≤ env (fun z => b z + y) x := by
    apply le_on_Icc_of_Ioo (f := fun x => b x - y) (hb.sub continuousOn_const) hcont
    intro x₀ hx₀
    apply le_env_of_chords hx₀ (by show b x₀ - y ≤ b x₀ + y; linarith)
    intro s x' hs hsx hxx' hx'
    exact chord_cond_of_gap hb hG hs hsx hxx' hx'
  intro x hx
  rw [abs_le]
  constructor
  · linarith [hlow x hx]
  · linarith [env_le hm hx]

/-- The envelope of a continuous `h` takes the value `h 1` at the endpoint 1 (cone argument). -/
theorem env_one {h : ℝ → ℝ} (hh : ContinuousOn h (Icc 0 1)) {m : ℝ}
    (hm : ∀ y ∈ Icc (0 : ℝ) 1, m ≤ h y) : env h 1 = h 1 := by
  have h1I : (1 : ℝ) ∈ Icc (0 : ℝ) 1 := ⟨zero_le_one, le_rfl⟩
  apply le_antisymm (env_le hm h1I)
  by_contra hcon
  push Not at hcon
  set η := (h 1 - env h 1) / 2 with hη
  have hη0 : 0 < η := by rw [hη]; linarith
  obtain ⟨ρ, hρ, hρh⟩ := Metric.continuousWithinAt_iff.mp (hh 1 h1I) η hη0
  set K := max 0 ((h 1 - η - m) / ρ) with hK
  have hK0 : 0 ≤ K := le_max_left _ _
  have hp : (K, h 1 - η - K) ∈ AffMin h := by
    intro y hy
    simp only
    by_cases hc : 1 - y < ρ
    · have hd : dist y 1 < ρ := by
        rw [Real.dist_eq, abs_sub_comm, abs_of_nonneg (by linarith [hy.2])]; exact hc
      have := hρh hy hd
      rw [Real.dist_eq, abs_sub_lt_iff] at this
      have : K * (y - 1) ≤ 0 := mul_nonpos_of_nonneg_of_nonpos hK0 (by linarith [hy.2])
      nlinarith
    · push Not at hc
      have hKρ : (h 1 - η - m) / ρ ≤ K := le_max_right _ _
      rw [div_le_iff₀ hρ] at hKρ
      have : K * ρ ≤ K * (1 - y) := mul_le_mul_of_nonneg_left hc hK0
      have := hm y hy
      nlinarith
  have := le_env h1I hp
  simp only at this
  linarith

/-- Endpoint terms rewritten (DT10 Lemma Λ0(i) form). -/
lemma endTerm_eq (b : ℝ → ℝ) (y : ℝ) {s x : ℝ} (hsx : s < x) :
    endTerm b y s x = b x + (1 - x) * (b x - b s - 2 * y) / (x - s) - y := by
  have hd : x - s ≠ 0 := by linarith
  unfold endTerm; field_simp; ring

lemma endSet_nonempty (b : ℝ → ℝ) (y : ℝ) : (endSet b y).Nonempty :=
  ⟨_, 0, 1 / 2, le_rfl, by norm_num, by norm_num, rfl⟩

/-- The endpoint formula is finite for `y > 0` (uniform continuity of `b`). -/
lemma endSet_bdd {b : ℝ → ℝ} (hb : ContinuousOn b (Icc 0 1)) {y : ℝ} (hy : 0 < y) :
    BddAbove (endSet b y) := by
  obtain ⟨C, hC⟩ := bound_of_cont hb
  have huc := isCompact_Icc.uniformContinuousOn_of_continuous hb
  obtain ⟨δ, hδ, hδb⟩ := Metric.uniformContinuousOn_iff.mp huc y hy
  refine ⟨C + 2 * (C + y) / δ, ?_⟩
  rintro v ⟨s, x, hs, hsx, hx1, rfl⟩
  rw [endTerm_eq b y hsx]
  have hsI : s ∈ Icc (0 : ℝ) 1 := ⟨hs, by linarith⟩
  have hxI : x ∈ Icc (0 : ℝ) 1 := ⟨by linarith, hx1.le⟩
  have hbx := abs_le.mp (hC x hxI)
  have hbs := abs_le.mp (hC s hsI)
  have hd : 0 < x - s := by linarith
  have h1x : 0 ≤ 1 - x := by linarith
  have hC0 : 0 ≤ C := (abs_nonneg _).trans (hC x hxI)
  by_cases hc : x - s < δ
  · have hdist : dist x s < δ := by rw [Real.dist_eq, abs_of_pos hd]; exact hc
    have := hδb x hxI s hsI hdist
    rw [Real.dist_eq, abs_sub_lt_iff] at this
    have hneg : (1 - x) * (b x - b s - 2 * y) / (x - s) ≤ 0 :=
      div_nonpos_of_nonpos_of_nonneg (mul_nonpos_of_nonneg_of_nonpos h1x (by linarith)) hd.le
    have : 0 ≤ 2 * (C + y) / δ := by positivity
    linarith
  · push Not at hc
    have hnum : (1 - x) * (b x - b s - 2 * y) ≤ 2 * (C + y) := by
      by_cases hn : b x - b s - 2 * y ≤ 0
      · have := mul_nonpos_of_nonneg_of_nonpos h1x hn; nlinarith
      · push Not at hn
        have : (1 - x) * (b x - b s - 2 * y) ≤ 1 * (b x - b s - 2 * y) :=
          mul_le_mul_of_nonneg_right (by linarith) hn.le
        linarith
    have hfrac : (1 - x) * (b x - b s - 2 * y) / (x - s) ≤ 2 * (C + y) / δ := by
      rcases le_or_gt ((1 - x) * (b x - b s - 2 * y)) 0 with hn | hn
      · have : (1 - x) * (b x - b s - 2 * y) / (x - s) ≤ 0 := div_nonpos_of_nonpos_of_nonneg hn hd.le
        have : 0 ≤ 2 * (C + y) / δ := by positivity
        linarith
      · calc (1 - x) * (b x - b s - 2 * y) / (x - s) ≤ (1 - x) * (b x - b s - 2 * y) / δ :=
              div_le_div_of_nonneg_left hn.le hδ hc
          _ ≤ 2 * (C + y) / δ := div_le_div_of_nonneg_right hnum hδ.le
    linarith

lemma endTerm_le_Lhat {b : ℝ → ℝ} (hb : ContinuousOn b (Icc 0 1)) {y : ℝ} (hy : 0 < y)
    {s x : ℝ} (hs : 0 ≤ s) (hsx : s < x) (hx1 : x < 1) : endTerm b y s x ≤ Lhat b y :=
  (le_csSup (endSet_bdd hb hy) ⟨s, x, hs, hsx, hx1, rfl⟩).trans (le_max_right _ _)

/-- Lemma F(c), necessity: every `y`-fit has endpoint `φ 1 ≥ L̂_b(y)`. -/
theorem Lhat_le_fit {b φ : ℝ → ℝ} {y : ℝ} (hf : IsFit b y φ) : Lhat b y ≤ φ 1 := by
  obtain ⟨_, hcv, hab⟩ := hf
  apply max_le
  · have := abs_le.mp (hab 1 ⟨zero_le_one, le_rfl⟩); linarith
  apply csSup_le (endSet_nonempty b y)
  rintro v ⟨s, x, hs, hsx, hx1, rfl⟩
  have hc := convex_le_chord hcv hs hsx hx1 le_rfl
  have h1 := abs_le.mp (hab x ⟨by linarith, hx1.le⟩)
  have h2 := abs_le.mp (hab s ⟨hs, by linarith⟩)
  have hd : 0 < x - s := by linarith
  have hd1 : 0 < 1 - s := by linarith
  unfold chord at hc
  unfold endTerm
  rw [div_le_iff₀ hd]
  have e : (1 - x) / (1 - s) * φ s + (x - s) / (1 - s) * φ 1
      = ((1 - x) * φ s + (x - s) * φ 1) / (1 - s) := by field_simp
  rw [e, le_div_iff₀ hd1] at hc
  have h1x : 0 ≤ 1 - x := by linarith
  nlinarith [mul_le_mul_of_nonneg_left h2.2 h1x]

/-- The maximal fit `env(b + y)` has endpoint value `b 1 + y`. -/
theorem fit_env_one {b : ℝ → ℝ} (hb : ContinuousOn b (Icc 0 1)) (y : ℝ) :
    env (fun z => b z + y) 1 = b 1 + y := by
  obtain ⟨C, hC⟩ := bound_of_cont hb
  have hm : ∀ z ∈ Icc (0 : ℝ) 1, -C + y ≤ b z + y := fun z hz => by
    have := neg_abs_le (b z); have := hC z hz; linarith
  exact env_one (hb.add continuousOn_const) hm

lemma Lhat_le_top {b : ℝ → ℝ} (hb : ContinuousOn b (Icc 0 1)) {y : ℝ} (hy : 0 ≤ y)
    (hG : Gch b ≤ 2 * y) : Lhat b y ≤ b 1 + y :=
  (Lhat_le_fit (fit_env hb hy hG)).trans (fit_env_one hb y).le

/-- Lemma F(c), sufficiency: a `y`-fit with endpoint exactly `L̂_b(y)` exists (`y > 0`). -/
theorem exists_fit_endpoint_eq {b : ℝ → ℝ} (hb : ContinuousOn b (Icc 0 1)) {y : ℝ} (hy : 0 < y)
    (hG : Gch b ≤ 2 * y) : ∃ φ, IsFit b y φ ∧ φ 1 = Lhat b y := by
  classical
  obtain ⟨C, hC⟩ := bound_of_cont hb
  set z := Lhat b y with hz
  have hztop : z ≤ b 1 + y := Lhat_le_top hb hy.le hG
  set h : ℝ → ℝ := fun x => if x = 1 then z else b x + y with hh
  have hzlow : b 1 - y ≤ z := le_max_left _ _
  have hm : ∀ w ∈ Icc (0 : ℝ) 1, min z (-C + y) ≤ h w := by
    intro w hw
    simp only [hh]
    split_ifs
    · exact min_le_left _ _
    · have := neg_abs_le (b w); have := hC w hw
      exact (min_le_right _ _).trans (by linarith)
  have hcont := env_continuousOn hm
  have hup : ∀ x ∈ Icc (0 : ℝ) 1, env h x ≤ b x + y := by
    intro x hx
    have := env_le hm hx
    simp only [hh] at this
    split_ifs at this with h1
    · subst h1; linarith
    · exact this
  have hlow : ∀ x ∈ Icc (0 : ℝ) 1, b x - y ≤ env h x := by
    apply le_on_Icc_of_Ioo (f := fun x => b x - y) (hb.sub continuousOn_const) hcont
    intro x₀ hx₀
    have hx₀1 : x₀ ≠ 1 := ne_of_lt hx₀.2
    have hhx₀ : h x₀ = b x₀ + y := by simp [hh, hx₀1]
    apply le_env_of_chords hx₀ (by show b x₀ - y ≤ h x₀; rw [hhx₀]; linarith)
    intro s x' hs hsx hxx' hx'
    have hs1 : s ≠ 1 := by intro h; linarith
    have hhs : h s = b s + y := by simp [hh, hs1]
    show (b x₀ - y) * (x' - s) ≤ (x' - x₀) * h s + (x₀ - s) * h x'
    rw [hhs]
    rcases lt_or_eq_of_le hx' with hx'1 | hx'1
    · have hx'ne : x' ≠ 1 := ne_of_lt hx'1
      have hhx' : h x' = b x' + y := by simp [hh, hx'ne]
      rw [hhx']
      exact chord_cond_of_gap hb hG hs hsx hxx' hx'
    · subst hx'1
      have hh1 : h 1 = z := by simp [hh]
      rw [hh1]
      have he := endTerm_le_Lhat hb hy hs hsx hx₀.2
      rw [← hz] at he
      unfold endTerm at he
      have hdx : 0 < x₀ - s := by linarith
      rw [div_le_iff₀ hdx] at he
      nlinarith
  refine ⟨env h, ⟨hcont, env_convex hm, fun x hx => abs_le.mpr ⟨by linarith [hlow x hx],
    by linarith [hup x hx]⟩⟩, le_antisymm ?_ ?_⟩
  · have := env_le hm (show (1 : ℝ) ∈ Icc (0 : ℝ) 1 from ⟨zero_le_one, le_rfl⟩)
    simpa [hh] using this
  · exact Lhat_le_fit ⟨hcont, env_convex hm, fun x hx => abs_le.mpr ⟨by linarith [hlow x hx],
      by linarith [hup x hx]⟩⟩

/-- Convex combinations of fits are fits. -/
lemma fit_combo {b φ ψ : ℝ → ℝ} {y lam : ℝ} (hφ : IsFit b y φ) (hψ : IsFit b y ψ)
    (h0 : 0 ≤ lam) (h1 : lam ≤ 1) : IsFit b y (fun x => lam * φ x + (1 - lam) * ψ x) := by
  obtain ⟨c1, v1, a1⟩ := hφ
  obtain ⟨c2, v2, a2⟩ := hψ
  refine ⟨(continuousOn_const.mul c1).add (continuousOn_const.mul c2),
    (v1.smul h0).add (v2.smul (by linarith)), fun x hx => ?_⟩
  have e1 := abs_le.mp (a1 x hx)
  have e2 := abs_le.mp (a2 x hx)
  rw [abs_le]
  constructor <;> nlinarith [mul_le_mul_of_nonneg_left e1.1 h0, mul_le_mul_of_nonneg_left e1.2 h0,
    mul_le_mul_of_nonneg_left e2.1 (by linarith : (0 : ℝ) ≤ 1 - lam),
    mul_le_mul_of_nonneg_left e2.2 (by linarith : (0 : ℝ) ≤ 1 - lam)]

/-- **Endpoint spectrum** (DT10 Lemma F(b)+(c)): for `y > 0` and `G(b) ≤ 2y`, the endpoint values
of `y`-fits are exactly `[L̂_b(y), b 1 + y]`. -/
theorem endpoint_spectrum {b : ℝ → ℝ} (hb : ContinuousOn b (Icc 0 1)) {y : ℝ} (hy : 0 < y)
    (hG : Gch b ≤ 2 * y) (z : ℝ) :
    (∃ φ, IsFit b y φ ∧ φ 1 = z) ↔ z ∈ Icc (Lhat b y) (b 1 + y) := by
  constructor
  · rintro ⟨φ, hf, rfl⟩
    refine ⟨Lhat_le_fit hf, ?_⟩
    have := abs_le.mp (hf.2.2 1 ⟨zero_le_one, le_rfl⟩); linarith
  · rintro ⟨hz1, hz2⟩
    obtain ⟨φ₀, hf₀, he₀⟩ := exists_fit_endpoint_eq hb hy hG
    have hF := fit_env hb hy.le hG
    have hF1 := fit_env_one hb y
    rcases eq_or_lt_of_le (hz1.trans hz2) with heq | hlt
    · exact ⟨φ₀, hf₀, by linarith⟩
    · set lam := (b 1 + y - z) / (b 1 + y - Lhat b y) with hlam
      have hd : 0 < b 1 + y - Lhat b y := by linarith
      have hl0 : 0 ≤ lam := div_nonneg (by linarith) hd.le
      have hl1 : lam ≤ 1 := by rw [div_le_one hd]; linarith
      refine ⟨_, fit_combo hf₀ hF hl0 hl1, ?_⟩
      simp only [he₀, hF1, hlam]
      field_simp
      ring

/-- Lemma F(a): a `y`-fit exists iff `G(b) ≤ 2y` (`y ≥ 0`). -/
theorem exists_fit_iff {b : ℝ → ℝ} (hb : ContinuousOn b (Icc 0 1)) {y : ℝ} (hy : 0 ≤ y) :
    (∃ φ, IsFit b y φ) ↔ Gch b ≤ 2 * y :=
  ⟨fun ⟨_, hf⟩ => Gch_le_of_fit hf, fun hG => ⟨_, fit_env hb hy hG⟩⟩

end Seam

open scoped Pointwise

namespace Seam

/-! ## Module A3: GCM identity (BP11 adapter), Λ0 (ii)-(iv), scaling, dependence on y -/

lemma lower_of_cont {b : ℝ → ℝ} (hb : ContinuousOn b (Icc 0 1)) :
    ∃ m, ∀ z ∈ Icc (0 : ℝ) 1, m ≤ b z := by
  obtain ⟨C, hC⟩ := bound_of_cont hb
  exact ⟨-C, fun z hz => by have := neg_abs_le (b z); have := hC z hz; linarith⟩

/-- `env b` is the GREATEST continuous convex minorant of `b` on `[0,1]` (the GCM of BP11, up to
the sign convention LCM(−b) = −GCM(b)). -/
theorem env_greatest {b ψ : ℝ → ℝ} (hb : ContinuousOn b (Icc 0 1))
    (hψc : ContinuousOn ψ (Icc 0 1)) (hψv : ConvexOn ℝ (Icc 0 1) ψ)
    (hψb : ∀ x ∈ Icc (0 : ℝ) 1, ψ x ≤ b x) : ∀ x ∈ Icc (0 : ℝ) 1, ψ x ≤ env b x := by
  obtain ⟨m, hm⟩ := lower_of_cont hb
  apply le_on_Icc_of_Ioo hψc (env_continuousOn hm)
  intro x₀ hx₀
  apply le_env_of_chords hx₀ (hψb x₀ ⟨hx₀.1.le, hx₀.2.le⟩)
  intro s x' hs hsx hxx' hx'
  have h1 := convex_le_chord hψv hs hsx hxx' hx'
  have h2 := chord_mono hsx hxx' (hψb s ⟨hs, by linarith⟩) (hψb x' ⟨by linarith, hx'⟩)
  have hd : 0 < x' - s := by linarith
  have h3 : ψ x₀ ≤ chord b s x' x₀ := h1.trans h2
  unfold chord at h3
  have e : (x' - x₀) / (x' - s) * b s + (x₀ - s) / (x' - s) * b x'
      = ((x' - x₀) * b s + (x₀ - s) * b x') / (x' - s) := by field_simp
  rw [e, le_div_iff₀ hd] at h3
  exact h3

/-- **Chord gap = envelope gap**: `Gch b = sup_{x∈[0,1]} (b x − env b x)` (least upper bound).
With `env_greatest`, this identifies `Gch` with the functional `sup(b − GCM b)` of BP11. -/
theorem Gch_isLUB_envGap {b : ℝ → ℝ} (hb : ContinuousOn b (Icc 0 1)) :
    IsLUB ((fun x => b x - env b x) '' Icc 0 1) (Gch b) := by
  obtain ⟨m, hm⟩ := lower_of_cont hb
  constructor
  · rintro _ ⟨x, hx, rfl⟩
    set g := Gch b / 2 with hg
    have hg0 : 0 ≤ g := by have := Gch_nonneg b; rw [hg]; linarith
    have hfit := fit_env hb hg0 (by rw [hg]; linarith)
    have hψ : ∀ z ∈ Icc (0 : ℝ) 1, env (fun z => b z + g) z - g ≤ b z := fun z hz => by
      have := abs_le.mp (hfit.2.2 z hz); linarith
    have hcv : ConvexOn ℝ (Icc 0 1) (fun z => env (fun z => b z + g) z - g) :=
      (hfit.2.1.add_const (-g)).congr (fun z _ => by simp [sub_eq_add_neg])
    have h1 : env (fun z => b z + g) x - g ≤ env b x :=
      env_greatest hb (hfit.1.sub continuousOn_const) hcv hψ x hx
    have h2 := abs_le.mp (hfit.2.2 x hx)
    show b x - env b x ≤ Gch b
    linarith
  · intro v hv
    have hI0 : (0 : ℝ) ∈ Icc (0 : ℝ) 1 := ⟨le_rfl, zero_le_one⟩
    apply max_le
    · have h : b 0 - env b 0 ≤ v := hv (mem_image_of_mem _ hI0)
      have := env_le hm hI0; linarith
    · apply csSup_le (gapSet_nonempty b)
      rintro _ ⟨s, x, x', hs, hsx, hxx', hx', rfl⟩
      have hxI : x ∈ Icc (0 : ℝ) 1 := ⟨by linarith, by linarith⟩
      have h1 := convex_le_chord (env_convex hm) hs hsx hxx' hx'
      have h2 := chord_mono hsx hxx' (env_le hm ⟨hs, by linarith⟩) (env_le hm ⟨by linarith, hx'⟩)
      have h : b x - env b x ≤ v := hv (mem_image_of_mem _ hxI)
      linarith

/-- DT11 Lemma 3(a) for `Gch`: the chord gap is at most the range. -/
theorem Gch_le_range {b : ℝ → ℝ} {w : ℝ}
    (hw : ∀ s ∈ Icc (0 : ℝ) 1, ∀ x ∈ Icc (0 : ℝ) 1, b x - b s ≤ w) : Gch b ≤ w := by
  have hI0 : (0 : ℝ) ∈ Icc (0 : ℝ) 1 := ⟨le_rfl, zero_le_one⟩
  have hw0 : 0 ≤ w := by have := hw 0 hI0 0 hI0; linarith
  apply max_le hw0
  apply csSup_le (gapSet_nonempty b)
  rintro _ ⟨s, x, x', hs, hsx, hxx', hx', rfl⟩
  obtain ⟨ha, hb', hab⟩ := chord_weights hsx hxx'
  have hxI : x ∈ Icc (0 : ℝ) 1 := ⟨by linarith, by linarith⟩
  have h1 := hw s ⟨hs, by linarith⟩ x hxI
  have h2 := hw x' ⟨by linarith, hx'⟩ x hxI
  unfold chord
  have k1 := mul_le_mul_of_nonneg_left h1 ha
  have k2 := mul_le_mul_of_nonneg_left h2 hb'
  have ex : b x = (x' - x) / (x' - s) * b x + (x - s) / (x' - s) * b x := by
    rw [← add_mul, hab, one_mul]
  have ew : w = (x' - x) / (x' - s) * w + (x - s) / (x' - s) * w := by
    rw [← add_mul, hab, one_mul]
  nlinarith

lemma gapSet_smul (b : ℝ → ℝ) (c : ℝ) : gapSet (fun x => c * b x) = c • gapSet b := by
  ext v
  rw [Set.mem_smul_set]
  constructor
  · rintro ⟨s, x, x', hs, hsx, hxx', hx', rfl⟩
    exact ⟨b x - chord b s x' x, ⟨s, x, x', hs, hsx, hxx', hx', rfl⟩,
      by simp only [smul_eq_mul, chord]; ring⟩
  · rintro ⟨w, ⟨s, x, x', hs, hsx, hxx', hx', rfl⟩, rfl⟩
    exact ⟨s, x, x', hs, hsx, hxx', hx', by simp only [smul_eq_mul, chord]; ring⟩

/-- Scaling: `G(c·b) = c·G(b)` for `c > 0`. -/
theorem Gch_smul (b : ℝ → ℝ) {c : ℝ} (hc : 0 < c) : Gch (fun x => c * b x) = c * Gch b := by
  unfold Gch
  rw [gapSet_smul, Real.sSup_smul_of_nonneg hc.le, smul_eq_mul, mul_max_of_nonneg _ _ hc.le,
    mul_zero]

lemma endSet_smul (b : ℝ → ℝ) (c y : ℝ) :
    endSet (fun x => c * b x) (c * y) = c • endSet b y := by
  ext v
  rw [Set.mem_smul_set]
  constructor
  · rintro ⟨s, x, hs, hsx, hx1, rfl⟩
    exact ⟨endTerm b y s x, ⟨s, x, hs, hsx, hx1, rfl⟩, by simp only [smul_eq_mul, endTerm]; ring⟩
  · rintro ⟨w, ⟨s, x, hs, hsx, hx1, rfl⟩, rfl⟩
    exact ⟨s, x, hs, hsx, hx1, by simp only [smul_eq_mul, endTerm]; ring⟩

/-- Scaling: `L̂_{c·b}(c·y) = c·L̂_b(y)` and `Λ_{c·b}(c·y) = c·Λ_b(y)` for `c > 0`. -/
theorem Lhat_smul (b : ℝ → ℝ) {c : ℝ} (hc : 0 < c) (y : ℝ) :
    Lhat (fun x => c * b x) (c * y) = c * Lhat b y := by
  unfold Lhat
  rw [endSet_smul, Real.sSup_smul_of_nonneg hc.le, smul_eq_mul, mul_max_of_nonneg _ _ hc.le]
  congr 1; ring

theorem Lam_smul (b : ℝ → ℝ) {c : ℝ} (hc : 0 < c) (y : ℝ) :
    Lam (fun x => c * b x) (c * y) = c * Lam b y := by
  unfold Lam; rw [Lhat_smul b hc]; ring

lemma Lhat_ge (b : ℝ → ℝ) (y : ℝ) : b 1 - y ≤ Lhat b y := le_max_left _ _

lemma Lam_ge (b : ℝ → ℝ) (y : ℝ) : b 1 ≤ Lam b y := by
  have := Lhat_ge b y; unfold Lam; linarith

/-- DT10 Lemma Λ0(ii): `y ↦ Λ_b(y)` is non-increasing (`y > 0`). -/
theorem Lam_antitone {b : ℝ → ℝ} (hb : ContinuousOn b (Icc 0 1)) {y y' : ℝ} (hy : 0 < y)
    (hyy : y ≤ y') : Lam b y' ≤ Lam b y := by
  unfold Lam Lhat
  have hy' : 0 < y' := by linarith
  have key : sSup (endSet b y') ≤ max (b 1 - y) (sSup (endSet b y)) + y - y' := by
    apply csSup_le (endSet_nonempty b y')
    rintro _ ⟨s, x, hs, hsx, hx1, rfl⟩
    have hle := endTerm_le_Lhat hb hy hs hsx hx1
    unfold Lhat at hle
    rw [endTerm_eq b y' hsx]
    rw [endTerm_eq b y hsx] at hle
    have hd : 0 < x - s := by linarith
    have h1x : 0 ≤ 1 - x := by linarith
    have : (1 - x) * (b x - b s - 2 * y') / (x - s) ≤ (1 - x) * (b x - b s - 2 * y) / (x - s) :=
      div_le_div_of_nonneg_right (mul_le_mul_of_nonneg_left (by linarith) h1x) hd.le
    linarith
  have h1 : b 1 - y' ≤ max (b 1 - y) (sSup (endSet b y)) + y - y' := by
    have := le_max_left (b 1 - y) (sSup (endSet b y)); linarith
  have := max_le h1 key
  linarith

/-- DT10 Lemma Λ0(iii): for `y ≥ G(b)/2`, `L̂_b(y) ≤ b 1 + G(b) − y` (so `Λ_b(y) ≤ b 1 + G(b)`). -/
theorem Lhat_le_G {b : ℝ → ℝ} (hb : ContinuousOn b (Icc 0 1)) {y : ℝ}
    (hG : Gch b ≤ 2 * y) : Lhat b y ≤ b 1 + Gch b - y := by
  set g := Gch b / 2 with hg
  have hg0 : 0 ≤ g := by have := Gch_nonneg b; rw [hg]; linarith
  have hfit := fit_env hb hg0 (by rw [hg]; linarith)
  have hφ : IsFit b y (fun z => env (fun z => b z + g) z + (g - y)) := by
    refine ⟨hfit.1.add continuousOn_const, hfit.2.1.add_const _, fun x hx => ?_⟩
    have := abs_le.mp (hfit.2.2 x hx)
    rw [abs_le]; constructor <;> linarith
  have h1 : Lhat b y ≤ env (fun z => b z + g) 1 + (g - y) := Lhat_le_fit hφ
  rw [fit_env_one hb g] at h1
  linarith

/-- DT10 Lemma Λ0(iv) (radius-shifted perturbation). -/
theorem Lhat_perturb {b b' : ℝ → ℝ} (hb' : ContinuousOn b' (Icc 0 1)) {y E : ℝ}
    (hE : ∀ x ∈ Icc (0 : ℝ) 1, |b x - b' x| ≤ E) (hyE : 0 < y - E)
    (hG : Gch b' ≤ 2 * (y - E)) : Lhat b y ≤ Lhat b' (y - E) := by
  obtain ⟨φ, hf, he⟩ := exists_fit_endpoint_eq hb' hyE hG
  have hfb : IsFit b y φ := by
    refine ⟨hf.1, hf.2.1, fun x hx => ?_⟩
    have h1 := abs_le.mp (hf.2.2 x hx)
    have h2 := abs_le.mp (hE x hx)
    rw [abs_le]; constructor <;> linarith
  rw [← he]; exact Lhat_le_fit hfb

theorem Lam_perturb {b b' : ℝ → ℝ} (hb' : ContinuousOn b' (Icc 0 1)) {y E : ℝ}
    (hE : ∀ x ∈ Icc (0 : ℝ) 1, |b x - b' x| ≤ E) (hyE : 0 < y - E)
    (hG : Gch b' ≤ 2 * (y - E)) : Lam b y ≤ Lam b' (y - E) + E := by
  have := Lhat_perturb hb' hE hyE hG; unfold Lam; linarith

/-- Perturbation of the chord gap: `G(b) ≤ G(b') + 2E` when `‖b − b'‖ ≤ E` on `[0,1]`. -/
theorem Gch_perturb {b b' : ℝ → ℝ} (hb' : ContinuousOn b' (Icc 0 1)) {E : ℝ}
    (hE : ∀ x ∈ Icc (0 : ℝ) 1, |b x - b' x| ≤ E) : Gch b ≤ Gch b' + 2 * E := by
  have hE0 : 0 ≤ E := (abs_nonneg _).trans (hE 0 ⟨le_rfl, zero_le_one⟩)
  apply max_le (by linarith [Gch_nonneg b'])
  apply csSup_le (gapSet_nonempty b)
  rintro _ ⟨s, x, x', hs, hsx, hxx', hx', rfl⟩
  have hg := gap_le_Gch hb' hs hsx hxx' hx'
  obtain ⟨ha, hb2, hab⟩ := chord_weights hsx hxx'
  have e1 := abs_le.mp (hE x ⟨by linarith, by linarith⟩)
  have e2 := abs_le.mp (hE s ⟨hs, by linarith⟩)
  have e3 := abs_le.mp (hE x' ⟨by linarith, hx'⟩)
  unfold chord at hg ⊢
  nlinarith [mul_le_mul_of_nonneg_left e2.1 ha, mul_le_mul_of_nonneg_left e3.1 hb2]

/-- `y ↦ L̂_b(y)` is convex on `(0, ∞)` (a sup of affine functions of `y`). -/
theorem Lhat_convex_y {b : ℝ → ℝ} (hb : ContinuousOn b (Icc 0 1)) {y₀ y₁ lam : ℝ}
    (h0 : 0 < y₀) (h1 : 0 < y₁) (hl0 : 0 ≤ lam) (hl1 : lam ≤ 1) :
    Lhat b ((1 - lam) * y₁ + lam * y₀) ≤ (1 - lam) * Lhat b y₁ + lam * Lhat b y₀ := by
  apply max_le
  · have a1 := Lhat_ge b y₁; have a0 := Lhat_ge b y₀
    nlinarith [mul_le_mul_of_nonneg_left a1 (by linarith : (0 : ℝ) ≤ 1 - lam),
      mul_le_mul_of_nonneg_left a0 hl0]
  · apply csSup_le (endSet_nonempty _ _)
    rintro _ ⟨s, x, hs, hsx, hx1, rfl⟩
    have e1 := endTerm_le_Lhat hb h1 hs hsx hx1
    have e0 := endTerm_le_Lhat hb h0 hs hsx hx1
    have hlin : endTerm b ((1 - lam) * y₁ + lam * y₀) s x
        = (1 - lam) * endTerm b y₁ s x + lam * endTerm b y₀ s x := by
      have hd : x - s ≠ 0 := by linarith
      unfold endTerm; field_simp; ring
    rw [hlin]
    exact add_le_add (mul_le_mul_of_nonneg_left e1 (by linarith))
      (mul_le_mul_of_nonneg_left e0 hl0)

/-- Right lower semicontinuity of `y ↦ L̂_b(y)` at `e > 0`, with a linear modulus. -/
theorem Lhat_right {b : ℝ → ℝ} (hb : ContinuousOn b (Icc 0 1)) {e γ : ℝ} (he : 0 < e)
    (hγ : 0 < γ) : ∃ K ≥ 0, ∀ δ ≥ 0, Lhat b e - γ - δ * K ≤ Lhat b (e + δ) := by
  have hlt : Lhat b e - γ < Lhat b e := by linarith
  unfold Lhat at hlt
  rcases lt_max_iff.mp hlt with h | h
  · refine ⟨1, zero_le_one, fun δ hδ => ?_⟩
    have := Lhat_ge b (e + δ)
    show max (b 1 - e) (sSup (endSet b e)) - γ - δ * 1 ≤ Lhat b (e + δ)
    linarith
  · obtain ⟨_, ⟨s, x, hs, hsx, hx1, rfl⟩, hv⟩ := exists_lt_of_lt_csSup (endSet_nonempty b e) h
    have hd : 0 < x - s := by linarith
    refine ⟨(2 - s - x) / (x - s), div_nonneg (by linarith) hd.le, fun δ hδ => ?_⟩
    have hle := endTerm_le_Lhat hb (by linarith : 0 < e + δ) hs hsx hx1
    have hlin : endTerm b (e + δ) s x = endTerm b e s x - δ * ((2 - s - x) / (x - s)) := by
      unfold endTerm; field_simp; ring
    rw [hlin] at hle
    show max (b 1 - e) (sSup (endSet b e)) - γ - δ * ((2 - s - x) / (x - s)) ≤ Lhat b (e + δ)
    linarith

/-- DT11 Lemma 3(b) via one affine fit `X − y + κ(1 − x)`: if `X` dominates `b(1−d) − κd`
(`hX`, with `d = 1 − x`), `X ≤ b x + w₀` (range ≤ `w₀`, `b 1 = 0`), `κ ≥ 0` and `w₀ + κ ≤ 2y`
(DT11: `κ = 2y − w₀`), then `Λ_b(y) ≤ X`. -/
theorem Lam_le_X {b : ℝ → ℝ} {y κ X w₀ : ℝ} (hκ : 0 ≤ κ) (hw : w₀ + κ ≤ 2 * y)
    (hX : ∀ x ∈ Icc (0 : ℝ) 1, b x - κ * (1 - x) ≤ X)
    (hXle : ∀ x ∈ Icc (0 : ℝ) 1, X ≤ b x + w₀) : Lam b y ≤ X := by
  have hfit : IsFit b y (fun x => X - y + κ * (1 - x)) := by
    refine ⟨by fun_prop, ⟨convex_Icc 0 1, fun x _ z _ a c ha hc hac => ?_⟩, fun x hx => ?_⟩
    · have : c = 1 - a := by linarith
      subst this
      simp only [smul_eq_mul]
      apply le_of_eq; ring
    · show |X - y + κ * (1 - x) - b x| ≤ y
      have h1 := hX x hx
      have h2 := hXle x hx
      have h3 : 0 ≤ κ * (1 - x) := mul_nonneg hκ (by linarith [hx.2])
      have h4 : κ * (1 - x) ≤ κ := by nlinarith [hx.1]
      rw [abs_le]; constructor <;> linarith
  have : Lhat b y ≤ X - y + κ * (1 - 1) := Lhat_le_fit hfit
  unfold Lam; linarith

end Seam

namespace Seam

/-! ## Module B: the seam (DT10 Lemma S(ii) + Lemma F(d)), distance, strict slack -/

/-- `K_u`: continuous on `[0,1]`, convex on `[0,u]`, concave on `[u,1]` (DT7-UNIF §0). -/
def Kset (u : ℝ) : Set (ℝ → ℝ) :=
  {φ | ContinuousOn φ (Icc 0 1) ∧ ConvexOn ℝ (Icc 0 u) φ ∧ ConcaveOn ℝ (Icc u 1) φ}

/-- A member of `K_u` within `e` of `B` on `[0,1]` exists. -/
def SeamFit (B : ℝ → ℝ) (u e : ℝ) : Prop :=
  ∃ φ ∈ Kset u, ∀ x ∈ Icc (0 : ℝ) 1, |B x - φ x| ≤ e

/-- `d_u(B) = d(B, K_u)` in the sup norm on `[0,1]`. -/
noncomputable def dK (u : ℝ) (B : ℝ → ℝ) : ℝ := sInf {e | 0 ≤ e ∧ SeamFit B u e}

/-- Left piece `a·b₁(s) = B(us) − s·B(u)` (DT10 Lemma S, `a = √u`). -/
def pieceL (B : ℝ → ℝ) (u : ℝ) (s : ℝ) : ℝ := B (u * s) - s * B u

/-- Right piece `ε·b₂(s) = s·B(u) − B(1 − (1−u)s)` (DT10 Lemma S/F(d), `ε = √(1−u)`). -/
def pieceR (B : ℝ → ℝ) (u : ℝ) (s : ℝ) : ℝ := s * B u - B (1 - (1 - u) * s)

lemma pieceL_one (B : ℝ → ℝ) (u : ℝ) : pieceL B u 1 = 0 := by simp [pieceL]

lemma pieceR_one (B : ℝ → ℝ) (u : ℝ) : pieceR B u 1 = 0 := by
  simp only [pieceR]; rw [show (1 : ℝ) - (1 - u) * 1 = u by ring]; ring

lemma pieceL_cont {B : ℝ → ℝ} (hB : ContinuousOn B (Icc 0 1)) {u : ℝ} (hu0 : 0 ≤ u)
    (hu1 : u ≤ 1) : ContinuousOn (pieceL B u) (Icc 0 1) := by
  unfold pieceL
  have hmap : MapsTo (fun s : ℝ => u * s) (Icc 0 1) (Icc 0 1) := fun s hs =>
    ⟨mul_nonneg hu0 hs.1, by nlinarith [hs.1, hs.2]⟩
  exact (hB.comp (by fun_prop) hmap).sub (by fun_prop)

lemma pieceR_cont {B : ℝ → ℝ} (hB : ContinuousOn B (Icc 0 1)) {u : ℝ} (hu0 : 0 ≤ u)
    (hu1 : u ≤ 1) : ContinuousOn (pieceR B u) (Icc 0 1) := by
  unfold pieceR
  have hmap : MapsTo (fun s : ℝ => 1 - (1 - u) * s) (Icc 0 1) (Icc 0 1) := fun s hs =>
    ⟨by nlinarith [hs.1, hs.2], by nlinarith [hs.1, hs.2]⟩
  exact (by fun_prop : ContinuousOn (fun s : ℝ => s * B u) (Icc 0 1)).sub (hB.comp (by fun_prop) hmap)

/-- Restriction (necessity half of the seam): a member of `K_u` within `e` of `B` gives `e`-fits
of both pieces whose endpoint values sum to 0. -/
theorem seam_restrict {B : ℝ → ℝ} {u e : ℝ} (hu0 : 0 < u) (hu1 : u < 1) (h : SeamFit B u e) :
    ∃ φ₁ φ₂, IsFit (pieceL B u) e φ₁ ∧ IsFit (pieceR B u) e φ₂ ∧ φ₁ 1 + φ₂ 1 = 0 := by
  obtain ⟨Φ, ⟨hc, hcv, hcc⟩, hfit⟩ := h
  refine ⟨fun s => Φ (u * s) - s * B u, fun s => -Φ (1 - (1 - u) * s) + s * B u, ⟨?_, ?_, ?_⟩,
    ⟨?_, ?_, ?_⟩, ?_⟩
  · have hmap : MapsTo (fun s : ℝ => u * s) (Icc 0 1) (Icc 0 1) := fun s hs =>
      ⟨mul_nonneg hu0.le hs.1, by nlinarith [hs.1, hs.2]⟩
    exact (hc.comp (by fun_prop) hmap).sub (by fun_prop)
  · refine ⟨convex_Icc 0 1, fun x hx z hz a c ha hc' hac => ?_⟩
    have hxI : u * x ∈ Icc (0 : ℝ) u := ⟨mul_nonneg hu0.le hx.1, by nlinarith [hx.1, hx.2]⟩
    have hzI : u * z ∈ Icc (0 : ℝ) u := ⟨mul_nonneg hu0.le hz.1, by nlinarith [hz.1, hz.2]⟩
    have h := hcv.2 hxI hzI ha hc' hac
    simp only [smul_eq_mul] at h ⊢
    have e1 : a * (u * x) + c * (u * z) = u * (a * x + c * z) := by ring
    rw [e1] at h
    nlinarith
  · intro s hs
    have hsI : u * s ∈ Icc (0 : ℝ) 1 := ⟨mul_nonneg hu0.le hs.1, by nlinarith [hs.1, hs.2]⟩
    have := hfit (u * s) hsI
    show |Φ (u * s) - s * B u - (B (u * s) - s * B u)| ≤ e
    rw [abs_sub_comm] at this
    convert this using 2; ring
  · have hmap : MapsTo (fun s : ℝ => 1 - (1 - u) * s) (Icc 0 1) (Icc 0 1) := fun s hs =>
      ⟨by nlinarith [hs.1, hs.2], by nlinarith [hs.1, hs.2]⟩
    exact (hc.comp (by fun_prop) hmap).neg.add (by fun_prop)
  · refine ⟨convex_Icc 0 1, fun x hx z hz a c ha hc' hac => ?_⟩
    have hxI : 1 - (1 - u) * x ∈ Icc u 1 := ⟨by nlinarith [hx.1, hx.2], by nlinarith [hx.1, hx.2]⟩
    have hzI : 1 - (1 - u) * z ∈ Icc u 1 := ⟨by nlinarith [hz.1, hz.2], by nlinarith [hz.1, hz.2]⟩
    have h := hcc.2 hxI hzI ha hc' hac
    simp only [smul_eq_mul] at h ⊢
    have e1 : a * (1 - (1 - u) * x) + c * (1 - (1 - u) * z) = 1 - (1 - u) * (a * x + c * z) := by
      linear_combination hac
    rw [e1] at h
    nlinarith
  · intro s hs
    have hsI : 1 - (1 - u) * s ∈ Icc (0 : ℝ) 1 := ⟨by nlinarith [hs.1, hs.2], by nlinarith [hs.1, hs.2]⟩
    have := hfit _ hsI
    show |-Φ (1 - (1 - u) * s) + s * B u - (s * B u - B (1 - (1 - u) * s))| ≤ e
    convert this using 2; ring
  · show Φ (u * 1) - 1 * B u + (-Φ (1 - (1 - u) * 1) + 1 * B u) = 0
    rw [show (1 : ℝ) - (1 - u) * 1 = u by ring, mul_one]; ring

/-- Gluing (sufficiency half; DT10 Lemma S(ii) + F(d)): `e`-fits of both pieces with endpoint
values summing to 0 give a member of `K_u` within `e` of `B`. -/
theorem seam_glue {B : ℝ → ℝ} {u e : ℝ} (hu0 : 0 < u) (hu1 : u < 1) {φ₁ φ₂ : ℝ → ℝ}
    (h₁ : IsFit (pieceL B u) e φ₁) (h₂ : IsFit (pieceR B u) e φ₂) (hsum : φ₁ 1 + φ₂ 1 = 0) :
    SeamFit B u e := by
  classical
  obtain ⟨c₁, v₁, f₁⟩ := h₁
  obtain ⟨c₂, v₂, f₂⟩ := h₂
  have hu1' : 0 < 1 - u := by linarith
  set ΦL : ℝ → ℝ := fun x => φ₁ (x / u) + x / u * B u with hΦL
  set ΦR : ℝ → ℝ := fun x => -φ₂ ((1 - x) / (1 - u)) + (1 - x) / (1 - u) * B u with hΦR
  set Φ : ℝ → ℝ := fun x => if x ≤ u then ΦL x else ΦR x with hΦ
  have hmatch : ΦL u = ΦR u := by
    simp only [hΦL, hΦR, div_self hu0.ne', div_self hu1'.ne']
    linarith
  have hEqL : EqOn Φ ΦL (Icc 0 u) := fun x hx => by simp [hΦ, hx.2]
  have hEqR : EqOn Φ ΦR (Icc u 1) := fun x hx => by
    simp only [hΦ]
    split_ifs with h
    · have : x = u := le_antisymm h hx.1
      rw [this]; exact hmatch
    · rfl
  have hmapL : MapsTo (fun x => x / u) (Icc 0 u) (Icc 0 1) := fun x hx =>
    ⟨div_nonneg hx.1 hu0.le, (div_le_one hu0).mpr hx.2⟩
  have hmapR : MapsTo (fun x => (1 - x) / (1 - u)) (Icc u 1) (Icc 0 1) := fun x hx =>
    ⟨div_nonneg (by linarith [hx.2]) hu1'.le, (div_le_one hu1').mpr (by linarith [hx.1])⟩
  have hcL : ContinuousOn ΦL (Icc 0 u) :=
    (c₁.comp (continuousOn_id.div_const u) hmapL).add
      ((continuousOn_id.div_const u).mul continuousOn_const)
  have hcR : ContinuousOn ΦR (Icc u 1) :=
    ((c₂.comp ((continuousOn_const.sub continuousOn_id).div_const (1 - u)) hmapR).neg).add
      (((continuousOn_const.sub continuousOn_id).div_const (1 - u)).mul continuousOn_const)
  refine ⟨Φ, ⟨?_, ?_, ?_⟩, ?_⟩
  · rw [← Icc_union_Icc_eq_Icc hu0.le hu1.le]
    exact (hcL.congr hEqL).union_of_isClosed (hcR.congr hEqR) isClosed_Icc isClosed_Icc
  · refine ConvexOn.congr ⟨convex_Icc 0 u, fun x hx z hz a c ha hc hac => ?_⟩ hEqL.symm
    have h := v₁.2 (hmapL hx) (hmapL hz) ha hc hac
    simp only [smul_eq_mul, hΦL] at h ⊢
    have e1 : (a * x + c * z) / u = a * (x / u) + c * (z / u) := by field_simp
    rw [e1]
    have e2 : (a * (x / u) + c * (z / u)) * B u = a * (x / u * B u) + c * (z / u * B u) := by ring
    rw [e2]
    nlinarith
  · refine ConcaveOn.congr ⟨convex_Icc u 1, fun x hx z hz a c ha hc hac => ?_⟩ hEqR.symm
    have h := v₂.2 (hmapR hx) (hmapR hz) ha hc hac
    simp only [smul_eq_mul, hΦR] at h ⊢
    have e1 : (1 - (a * x + c * z)) / (1 - u) = a * ((1 - x) / (1 - u)) + c * ((1 - z) / (1 - u)) := by
      have hc' : c = 1 - a := by linarith
      rw [hc']; field_simp; ring
    rw [e1]
    have e2 : (a * ((1 - x) / (1 - u)) + c * ((1 - z) / (1 - u))) * B u
        = a * ((1 - x) / (1 - u) * B u) + c * ((1 - z) / (1 - u) * B u) := by ring
    rw [e2]
    nlinarith
  · intro x hx
    by_cases hxu : x ≤ u
    · rw [hEqL ⟨hx.1, hxu⟩]
      have := f₁ (x / u) (hmapL ⟨hx.1, hxu⟩)
      simp only [hΦL, pieceL] at this ⊢
      rw [show u * (x / u) = x by field_simp] at this
      rw [abs_sub_comm]
      convert this using 2; ring
    · push Not at hxu
      rw [hEqR ⟨hxu.le, hx.2⟩]
      have := f₂ ((1 - x) / (1 - u)) (hmapR ⟨hxu.le, hx.2⟩)
      simp only [hΦR, pieceR] at this ⊢
      rw [show 1 - (1 - u) * ((1 - x) / (1 - u)) = x by field_simp; ring] at this
      convert this using 2; ring

/-- **Lemma F(d) at the feasibility level** (exact radius `e > 0`): a member of `K_u` within `e`
exists iff `G(Aₗ) ≤ 2e`, `G(Aᵣ) ≤ 2e` and `L̂_{Aₗ}(e) + L̂_{Aᵣ}(e) ≤ 0` (scaled pieces). -/
theorem seamFit_iff {B : ℝ → ℝ} (hB : ContinuousOn B (Icc 0 1)) {u e : ℝ} (hu0 : 0 < u)
    (hu1 : u < 1) (he : 0 < e) :
    SeamFit B u e ↔ Gch (pieceL B u) ≤ 2 * e ∧ Gch (pieceR B u) ≤ 2 * e ∧
      Lhat (pieceL B u) e + Lhat (pieceR B u) e ≤ 0 := by
  have hL := pieceL_cont hB hu0.le hu1.le
  have hR := pieceR_cont hB hu0.le hu1.le
  constructor
  · intro h
    obtain ⟨φ₁, φ₂, f₁, f₂, hs⟩ := seam_restrict hu0 hu1 h
    exact ⟨Gch_le_of_fit f₁, Gch_le_of_fit f₂, by
      have := Lhat_le_fit f₁; have := Lhat_le_fit f₂; linarith⟩
  · rintro ⟨hG₁, hG₂, hsum⟩
    have t₁ := Lhat_le_top hL he.le hG₁
    have t₂ := Lhat_le_top hR he.le hG₂
    rw [pieceL_one] at t₁; rw [pieceR_one] at t₂
    set v := max (Lhat (pieceL B u) e) (-e) with hv
    obtain ⟨φ₁, f₁, e₁⟩ := (endpoint_spectrum hL he hG₁ v).mpr
      ⟨le_max_left _ _, by rw [pieceL_one]; exact max_le (by linarith) (by linarith)⟩
    obtain ⟨φ₂, f₂, e₂⟩ := (endpoint_spectrum hR he hG₂ (-v)).mpr
      ⟨by rw [hv]; rcases le_total (Lhat (pieceL B u) e) (-e) with h | h
          · rw [max_eq_right h]; linarith
          · rw [max_eq_left h]; linarith,
       by rw [pieceR_one, hv]; have := le_max_right (Lhat (pieceL B u) e) (-e); linarith⟩
    exact seam_glue hu0 hu1 f₁ f₂ (by rw [e₁, e₂]; ring)

lemma seamFit_mono {B : ℝ → ℝ} {u e e' : ℝ} (h : SeamFit B u e) (hee : e ≤ e') :
    SeamFit B u e' := by
  obtain ⟨φ, hφ, hf⟩ := h
  exact ⟨φ, hφ, fun x hx => (hf x hx).trans hee⟩

lemma zero_mem_Kset (u : ℝ) : (fun _ => (0 : ℝ)) ∈ Kset u :=
  ⟨continuousOn_const, convexOn_const 0 (convex_Icc 0 u), concaveOn_const 0 (convex_Icc u 1)⟩

lemma seamFit_exists {B : ℝ → ℝ} (hB : ContinuousOn B (Icc 0 1)) (u : ℝ) :
    ∃ e, 0 ≤ e ∧ SeamFit B u e := by
  obtain ⟨C, hC⟩ := bound_of_cont hB
  have hC0 : 0 ≤ C := (abs_nonneg _).trans (hC 0 ⟨le_rfl, zero_le_one⟩)
  exact ⟨C, hC0, _, zero_mem_Kset u, fun x hx => by simpa using hC x hx⟩

lemma dK_set_bdd (B : ℝ → ℝ) (u : ℝ) : BddBelow {e | 0 ≤ e ∧ SeamFit B u e} :=
  ⟨0, fun _ h => h.1⟩

lemma dK_nonneg (B : ℝ → ℝ) (u : ℝ) : 0 ≤ dK u B := by
  unfold dK
  rcases Set.eq_empty_or_nonempty {e | 0 ≤ e ∧ SeamFit B u e} with h | h
  · rw [h, Real.sInf_empty]
  · exact le_csInf h (fun _ he => he.1)

/-- From `d_u(B) < e'` (or `≤ e`), fits at every larger radius. -/
lemma seamFit_of_dK_lt {B : ℝ → ℝ} (hB : ContinuousOn B (Icc 0 1)) {u e : ℝ}
    (h : dK u B < e) : SeamFit B u e := by
  have hne : ({e | 0 ≤ e ∧ SeamFit B u e} : Set ℝ).Nonempty := by
    obtain ⟨e, h0, hf⟩ := seamFit_exists hB u; exact ⟨e, h0, hf⟩
  obtain ⟨e₀, he₀, hlt⟩ := exists_lt_of_csInf_lt hne h
  exact seamFit_mono he₀.2 hlt.le

/-- **Distance = feasibility at the same radius** (named theorem):
for `0 < u < 1` and `e > 0`, `d_u(B) ≤ e ↔ SeamFit B u e` (the infimum is attained). -/
theorem dK_le_iff {B : ℝ → ℝ} (hB : ContinuousOn B (Icc 0 1)) {u e : ℝ} (hu0 : 0 < u)
    (hu1 : u < 1) (he : 0 < e) : dK u B ≤ e ↔ SeamFit B u e := by
  constructor
  · intro hd
    have hL := pieceL_cont hB hu0.le hu1.le
    have hR := pieceR_cont hB hu0.le hu1.le
    have hδ : ∀ δ > 0, Gch (pieceL B u) ≤ 2 * (e + δ) ∧ Gch (pieceR B u) ≤ 2 * (e + δ) ∧
        Lhat (pieceL B u) (e + δ) + Lhat (pieceR B u) (e + δ) ≤ 0 := fun δ hδ =>
      (seamFit_iff hB hu0 hu1 (by linarith)).mp (seamFit_of_dK_lt hB (by linarith))
    rw [seamFit_iff hB hu0 hu1 he]
    have hG : ∀ G, (∀ δ > 0, G ≤ 2 * (e + δ)) → G ≤ 2 * e := fun G hG => by
      by_contra hc; push Not at hc
      have := hG ((G - 2 * e) / 4) (by linarith); linarith
    refine ⟨hG _ (fun δ hδ' => (hδ δ hδ').1), hG _ (fun δ hδ' => (hδ δ hδ').2.1), ?_⟩
    by_contra hc; push Not at hc
    set σ := Lhat (pieceL B u) e + Lhat (pieceR B u) e with hσ
    obtain ⟨K₁, hK₁, h₁⟩ := Lhat_right hL he (show 0 < σ / 4 by linarith)
    obtain ⟨K₂, hK₂, h₂⟩ := Lhat_right hR he (show 0 < σ / 4 by linarith)
    set δ := σ / (4 * (K₁ + K₂ + 1)) with hδdef
    have hδ0 : 0 < δ := by positivity
    have a₁ := h₁ δ hδ0.le
    have a₂ := h₂ δ hδ0.le
    have a₃ := (hδ δ hδ0).2.2
    have hK : δ * K₁ + δ * K₂ ≤ σ / 4 := by
      have : δ * (K₁ + K₂ + 1) = σ / 4 := by rw [hδdef]; field_simp
      nlinarith
    linarith
  · intro h
    exact csInf_le (dK_set_bdd B u) ⟨he.le, h⟩

/-- **Lemma F(d)** (exact form, DT10 §2): for `0 < u < 1`, `e > 0`,
`d_u(B) ≤ e ↔ G(Aₗ) ≤ 2e ∧ G(Aᵣ) ≤ 2e ∧ L̂_{Aₗ}(e) + L̂_{Aᵣ}(e) ≤ 0`. -/
theorem lemmaF_d {B : ℝ → ℝ} (hB : ContinuousOn B (Icc 0 1)) {u e : ℝ} (hu0 : 0 < u)
    (hu1 : u < 1) (he : 0 < e) :
    dK u B ≤ e ↔ Gch (pieceL B u) ≤ 2 * e ∧ Gch (pieceR B u) ≤ 2 * e ∧
      Lhat (pieceL B u) e + Lhat (pieceR B u) e ≤ 0 :=
  (dK_le_iff hB hu0 hu1 he).trans (seamFit_iff hB hu0 hu1 he)

/-- **Strict slack**: `G(Aᵣ) < 2t` and `G(Aₗ) + Λ_{Aᵣ}(t) < 2t` give `d_u(B) < t`. -/
theorem strict_slack {B : ℝ → ℝ} (hB : ContinuousOn B (Icc 0 1)) {u t : ℝ} (hu0 : 0 < u)
    (hu1 : u < 1) (ht : 0 < t) (h₂ : Gch (pieceR B u) < 2 * t)
    (h₁ : Gch (pieceL B u) + Lam (pieceR B u) t < 2 * t) : dK u B < t := by
  have hL := pieceL_cont hB hu0.le hu1.le
  have hR := pieceR_cont hB hu0.le hu1.le
  set gL := Gch (pieceL B u)
  set gR := Gch (pieceR B u)
  set s₁ := 2 * t - (gL + Lam (pieceR B u) t) with hs₁
  set s₂ := 2 * t - gR with hs₂
  have hs₁0 : 0 < s₁ := by linarith
  have hs₂0 : 0 < s₂ := by linarith
  set D := |Lhat (pieceR B u) (t / 2) - Lhat (pieceR B u) t| with hD
  have hD0 : 0 ≤ D := abs_nonneg _
  set lam := min (1 / 2) (min (s₁ / (t + D + 1)) (s₂ / (t + 1))) with hlam
  have hl0 : 0 < lam := lt_min (by norm_num) (lt_min (by positivity) (by positivity))
  have hl1 : lam ≤ 1 / 2 := min_le_left _ _
  have hl2 : lam ≤ s₁ / (t + D + 1) := (min_le_right _ _).trans (min_le_left _ _)
  have hl3 : lam ≤ s₂ / (t + 1) := (min_le_right _ _).trans (min_le_right _ _)
  set t' := (1 - lam) * t + lam * (t / 2) with ht'
  have ht'lt : t' < t := by rw [ht']; nlinarith
  have ht'pos : 0 < t' := by rw [ht']; nlinarith
  have hlamt : lam * (t + D + 1) ≤ s₁ := by rwa [le_div_iff₀ (by positivity)] at hl2
  have hlamt2 : lam * (t + 1) ≤ s₂ := by rwa [le_div_iff₀ (by positivity)] at hl3
  have hGLnn := Gch_nonneg (pieceL B u)
  have hLamR : 0 ≤ Lam (pieceR B u) t := by have := Lam_ge (pieceR B u) t; rwa [pieceR_one] at this
  have hGR' : gR ≤ 2 * t' := by rw [ht']; nlinarith
  have hGL' : gL ≤ 2 * t' := by rw [ht']; nlinarith
  have hconv := Lhat_convex_y hR (show 0 < t / 2 by linarith) ht hl0.le (by linarith)
  rw [← ht'] at hconv
  have hDle : Lhat (pieceR B u) (t / 2) - Lhat (pieceR B u) t ≤ D := le_abs_self _
  have hLL := Lhat_le_G hL hGL'
  rw [pieceL_one] at hLL
  have hsum : Lhat (pieceL B u) t' + Lhat (pieceR B u) t' ≤ 0 := by
    have : Lhat (pieceR B u) t' ≤ Lhat (pieceR B u) t + lam * D := by nlinarith
    unfold Lam at hs₁
    rw [ht'] at hLL ⊢
    nlinarith
  have hfit := (seamFit_iff hB hu0 hu1 ht'pos).mpr ⟨hGL', hGR', hsum⟩
  exact lt_of_le_of_lt (csInf_le (dK_set_bdd B u) ⟨ht'pos.le, hfit⟩) ht'lt

/-! ### The edge `u = 1`: `d_1(B) = G(B)/2` -/

lemma concaveOn_Icc_one (φ : ℝ → ℝ) : ConcaveOn ℝ (Icc 1 1) φ := by
  refine ⟨convex_Icc 1 1, fun x hx z hz a c _ _ hac => ?_⟩
  have hx1 : x = 1 := le_antisymm hx.2 hx.1
  have hz1 : z = 1 := le_antisymm hz.2 hz.1
  subst hx1; subst hz1
  simp only [smul_eq_mul]
  have : a * 1 + c * 1 = 1 := by linarith
  rw [this, ← add_mul, hac, one_mul]

theorem seamFit_one_iff {B : ℝ → ℝ} (hB : ContinuousOn B (Icc 0 1)) {e : ℝ} (he : 0 ≤ e) :
    SeamFit B 1 e ↔ Gch B ≤ 2 * e := by
  rw [← exists_fit_iff hB he]
  constructor
  · rintro ⟨φ, ⟨hc, hv, _⟩, hf⟩
    exact ⟨φ, hc, hv, fun x hx => by rw [abs_sub_comm]; exact hf x hx⟩
  · rintro ⟨φ, hc, hv, hf⟩
    exact ⟨φ, ⟨hc, hv, concaveOn_Icc_one φ⟩, fun x hx => by rw [abs_sub_comm]; exact hf x hx⟩

/-- At the edge, `d_1(B) = G(B)/2` exactly (DT10 §0; BP11's `M/2`). -/
theorem dK_one {B : ℝ → ℝ} (hB : ContinuousOn B (Icc 0 1)) : dK 1 B = Gch B / 2 := by
  have hG0 := Gch_nonneg B
  apply le_antisymm
  · exact csInf_le (dK_set_bdd B 1) ⟨by linarith, (seamFit_one_iff hB (by linarith)).mpr (by linarith)⟩
  · apply le_csInf
    · exact ⟨_, by linarith, (seamFit_one_iff hB (by linarith)).mpr (le_refl (Gch B) |>.trans (by linarith))⟩
    · rintro e ⟨he0, hf⟩
      have := (seamFit_one_iff hB he0).mp hf; linarith

end Seam

namespace Seam

/-! ## Module C1: Lemma B, deterministic core, in the normalised bridges b₁, b₂ of DT10

 fix: `b₁`, `b₂` are the CLIPPED pieces (values on `[0,1]` only), so that Lemma S(i)
(independence of `b₁` and `b₂`) can be stated for them as `ℝ → ℝ`-valued maps. -/

/-- Clipping to `[0,1]`. -/
noncomputable def clip (s : ℝ) : ℝ := max 0 (min s 1)

lemma clip_of_mem {s : ℝ} (hs : s ∈ Icc (0 : ℝ) 1) : clip s = s := by
  unfold clip; rw [min_eq_left hs.2, max_eq_right hs.1]

lemma clip_mem (s : ℝ) : clip s ∈ Icc (0 : ℝ) 1 :=
  ⟨le_max_left _ _, max_le (by norm_num) (min_le_right _ _)⟩

lemma continuous_clip : Continuous clip := by unfold clip; fun_prop

lemma Gch_congr {f g : ℝ → ℝ} (h : ∀ x ∈ Icc (0 : ℝ) 1, f x = g x) : Gch f = Gch g := by
  unfold Gch; congr 2
  ext v
  constructor <;> rintro ⟨s, x, x', hs, hsx, hxx', hx', rfl⟩ <;>
    refine ⟨s, x, x', hs, hsx, hxx', hx', ?_⟩ <;>
    simp only [chord, h s ⟨hs, by linarith⟩, h x ⟨by linarith, by linarith⟩, h x' ⟨by linarith, hx'⟩]

lemma Lhat_congr {f g : ℝ → ℝ} (h : ∀ x ∈ Icc (0 : ℝ) 1, f x = g x) (y : ℝ) :
    Lhat f y = Lhat g y := by
  unfold Lhat
  rw [h 1 ⟨zero_le_one, le_rfl⟩]
  congr 2
  ext v
  constructor <;> rintro ⟨s, x, hs, hsx, hx1, rfl⟩ <;> refine ⟨s, x, hs, hsx, hx1, ?_⟩ <;>
    simp only [endTerm, h s ⟨hs, by linarith⟩, h x ⟨by linarith, hx1.le⟩]

lemma Lam_congr {f g : ℝ → ℝ} (h : ∀ x ∈ Icc (0 : ℝ) 1, f x = g x) (y : ℝ) :
    Lam f y = Lam g y := by unfold Lam; rw [Lhat_congr h]

/-- `b₁(s) = (B(us) − s·B(u))/√u` on `[0,1]` (DT10 Lemma S), clipped outside. -/
noncomputable def bL (B : ℝ → ℝ) (u : ℝ) (s : ℝ) : ℝ := pieceL B u (clip s) / √u

/-- `b₂(s) = −r₂(1 − s) = (s·B(u) − B(1 − (1−u)s))/√(1−u)` on `[0,1]` (DT10 Lemma F(d)),
clipped outside. -/
noncomputable def bR (B : ℝ → ℝ) (u : ℝ) (s : ℝ) : ℝ := pieceR B u (clip s) / √(1 - u)

lemma pieceL_eq (B : ℝ → ℝ) {u : ℝ} (hu0 : 0 < u) :
    ∀ s ∈ Icc (0 : ℝ) 1, pieceL B u s = √u * bL B u s := by
  intro s hs
  have : 0 < √u := Real.sqrt_pos.mpr hu0
  unfold bL; rw [clip_of_mem hs]; field_simp

lemma pieceR_eq (B : ℝ → ℝ) {u : ℝ} (hu1 : u < 1) :
    ∀ s ∈ Icc (0 : ℝ) 1, pieceR B u s = √(1 - u) * bR B u s := by
  intro s hs
  have : 0 < √(1 - u) := Real.sqrt_pos.mpr (by linarith)
  unfold bR; rw [clip_of_mem hs]; field_simp

lemma bL_cont {B : ℝ → ℝ} (hB : ContinuousOn B (Icc 0 1)) {u : ℝ} (hu0 : 0 < u) (hu1 : u < 1) :
    ContinuousOn (bL B u) (Icc 0 1) :=
  ((pieceL_cont hB hu0.le hu1.le).div_const (√u)).congr (fun s hs => by
    simp only [bL, clip_of_mem hs])

lemma bR_cont {B : ℝ → ℝ} (hB : ContinuousOn B (Icc 0 1)) {u : ℝ} (hu0 : 0 < u) (hu1 : u < 1) :
    ContinuousOn (bR B u) (Icc 0 1) :=
  ((pieceR_cont hB hu0.le hu1.le).div_const (√(1 - u))).congr (fun s hs => by
    simp only [bR, clip_of_mem hs])

lemma bR_one (B : ℝ → ℝ) (u : ℝ) : bR B u 1 = 0 := by
  simp [bR, clip_of_mem (show (1 : ℝ) ∈ Icc (0 : ℝ) 1 from ⟨zero_le_one, le_rfl⟩), pieceR_one]

lemma bL_one (B : ℝ → ℝ) (u : ℝ) : bL B u 1 = 0 := by
  simp [bL, clip_of_mem (show (1 : ℝ) ∈ Icc (0 : ℝ) 1 from ⟨zero_le_one, le_rfl⟩), pieceL_one]

lemma Gch_pieceL (B : ℝ → ℝ) {u : ℝ} (hu0 : 0 < u) :
    Gch (pieceL B u) = √u * Gch (bL B u) := by
  rw [Gch_congr (pieceL_eq B hu0), Gch_smul _ (Real.sqrt_pos.mpr hu0)]

lemma Gch_pieceR (B : ℝ → ℝ) {u : ℝ} (hu1 : u < 1) :
    Gch (pieceR B u) = √(1 - u) * Gch (bR B u) := by
  rw [Gch_congr (pieceR_eq B hu1), Gch_smul _ (Real.sqrt_pos.mpr (by linarith))]

lemma Lam_pieceR (B : ℝ → ℝ) {u : ℝ} (hu1 : u < 1) (t : ℝ) :
    Lam (pieceR B u) t = √(1 - u) * Lam (bR B u) (t / √(1 - u)) := by
  have hε : 0 < √(1 - u) := Real.sqrt_pos.mpr (by linarith)
  rw [Lam_congr (pieceR_eq B hu1), ← Lam_smul _ hε]; congr 1; field_simp

/-- **Lemma B, deterministic core** (DT10 §6, proved for the non-strict tail `d_u ≥ t`):
with `a = √u`, `ε = √(1−u)`, `G₁ = G(b₁)`, `G₂ = G(b₂)`, `Λ₂ = Λ_{b₂}(t/ε)`:
`t ≤ d_u(B)` implies `2t ≤ ε·G₂`, or `ε·G₂ < 2t` and `(2t − ε·Λ₂)/a ≤ G₁`. -/
theorem lemmaB_det {B : ℝ → ℝ} (hB : ContinuousOn B (Icc 0 1)) {u t : ℝ} (hu0 : 0 < u)
    (hu1 : u < 1) (ht : 0 < t) (hd : t ≤ dK u B) :
    2 * t ≤ √(1 - u) * Gch (bR B u) ∨
      (√(1 - u) * Gch (bR B u) < 2 * t ∧
        (2 * t - √(1 - u) * Lam (bR B u) (t / √(1 - u))) / √u ≤ Gch (bL B u)) := by
  have ha : 0 < √u := Real.sqrt_pos.mpr hu0
  by_contra hcon
  push Not at hcon
  obtain ⟨h₂, h₁⟩ := hcon
  have h₁' := h₁ h₂
  rw [lt_div_iff₀ ha] at h₁'
  have hs := strict_slack hB hu0 hu1 ht (by rw [Gch_pieceR B hu1]; exact h₂)
    (by rw [Gch_pieceL B hu0, Lam_pieceR B hu1]; nlinarith)
  linarith

end Seam

namespace Seam

/-! ## Module C2: rational-index surrogates, agreement on continuous paths, measurability
(the raw functionals are not product-measurable on all of ℝ → ℝ; the
surrogates are, and they agree with the raw functionals on paths continuous on [0,1].) -/

lemma le_of_near {α : Type*} [PseudoMetricSpace α] {G : α → ℝ} {p : α} (hG : ContinuousAt G p)
    {M : ℝ} (hM : ∀ ρ > 0, ∃ q, dist q p < ρ ∧ G q ≤ M) : G p ≤ M := by
  by_contra h
  push Not at h
  obtain ⟨ρ, hρ, hρG⟩ := Metric.continuousAt_iff.mp hG (G p - M) (by linarith)
  obtain ⟨q, hq, hqM⟩ := hM ρ hρ
  have := hρG hq
  rw [Real.dist_eq, abs_sub_lt_iff] at this
  linarith

lemma rat_triple {s x x' ρ : ℝ} (hs : 0 ≤ s) (hsx : s < x) (hxx' : x < x') (hx' : x' ≤ 1)
    (hρ : 0 < ρ) : ∃ q₁ q₂ q₃ : ℚ, 0 ≤ (q₁ : ℝ) ∧ (q₁ : ℝ) < q₂ ∧ (q₂ : ℝ) < q₃ ∧ (q₃ : ℝ) ≤ 1 ∧
      |(q₁ : ℝ) - s| < ρ ∧ |(q₂ : ℝ) - x| < ρ ∧ |(q₃ : ℝ) - x'| < ρ := by
  set r := min ρ (min ((x - s) / 3) ((x' - x) / 3)) with hr
  have hr0 : 0 < r := lt_min hρ (lt_min (by linarith) (by linarith))
  have hr1 : r ≤ ρ := min_le_left _ _
  have hr2 : r ≤ (x - s) / 3 := (min_le_right _ _).trans (min_le_left _ _)
  have hr3 : r ≤ (x' - x) / 3 := (min_le_right _ _).trans (min_le_right _ _)
  obtain ⟨q₁, h1a, h1b⟩ := exists_rat_btwn (show s < s + r by linarith)
  obtain ⟨q₂, h2a, h2b⟩ := exists_rat_btwn (show x < x + r by linarith)
  obtain ⟨q₃, h3a, h3b⟩ := exists_rat_btwn (show x' - r < x' by linarith)
  refine ⟨q₁, q₂, q₃, by linarith, by linarith, by linarith, by linarith, ?_, ?_, ?_⟩ <;>
    rw [abs_lt] <;> constructor <;> linarith

lemma rat_pair {s x ρ : ℝ} (hs : 0 ≤ s) (hsx : s < x) (hx1 : x < 1) (hρ : 0 < ρ) :
    ∃ q₁ q₂ : ℚ, 0 ≤ (q₁ : ℝ) ∧ (q₁ : ℝ) < q₂ ∧ (q₂ : ℝ) < 1 ∧
      |(q₁ : ℝ) - s| < ρ ∧ |(q₂ : ℝ) - x| < ρ := by
  set r := min ρ (min ((x - s) / 3) ((1 - x) / 2)) with hr
  have hr0 : 0 < r := lt_min hρ (lt_min (by linarith) (by linarith))
  have hr1 : r ≤ ρ := min_le_left _ _
  have hr2 : r ≤ (x - s) / 3 := (min_le_right _ _).trans (min_le_left _ _)
  have hr3 : r ≤ (1 - x) / 2 := (min_le_right _ _).trans (min_le_right _ _)
  obtain ⟨q₁, h1a, h1b⟩ := exists_rat_btwn (show s < s + r by linarith)
  obtain ⟨q₂, h2a, h2b⟩ := exists_rat_btwn (show x < x + r by linarith)
  refine ⟨q₁, q₂, by linarith, by linarith, by linarith, ?_, ?_⟩ <;>
    rw [abs_lt] <;> constructor <;> linarith

lemma rat_single {d ρ : ℝ} (hd : d ∈ Icc (0 : ℝ) 1) (hρ : 0 < ρ) :
    ∃ q : ℚ, (q : ℝ) ∈ Icc (0 : ℝ) 1 ∧ |(q : ℝ) - d| < ρ := by
  have h : max 0 (d - ρ) < min 1 (d + ρ) :=
    max_lt (lt_min (by norm_num) (by linarith [hd.1])) (lt_min (by linarith [hd.2]) (by linarith))
  obtain ⟨q, hqa, hqb⟩ := exists_rat_btwn h
  have h1 := le_max_left 0 (d - ρ)
  have h2 := le_max_right 0 (d - ρ)
  have h3 := min_le_left 1 (d + ρ)
  have h4 := min_le_right 1 (d + ρ)
  exact ⟨q, ⟨by linarith, by linarith⟩, by rw [abs_lt]; constructor <;> linarith⟩

/-! ### Index sets and surrogates -/

def Q3 : Set (ℚ × ℚ × ℚ) :=
  {q | 0 ≤ (q.1 : ℝ) ∧ (q.1 : ℝ) < q.2.1 ∧ (q.2.1 : ℝ) < q.2.2 ∧ (q.2.2 : ℝ) ≤ 1}
def Q2 : Set (ℚ × ℚ) := {q | 0 ≤ (q.1 : ℝ) ∧ (q.1 : ℝ) < q.2 ∧ (q.2 : ℝ) < 1}
def Q1 : Set ℚ := {q | 0 ≤ (q : ℝ) ∧ (q : ℝ) ≤ 1}

instance : Nonempty Q3 := ⟨⟨(0, 1 / 2, 1), by norm_num [Q3]⟩⟩
instance : Nonempty Q2 := ⟨⟨(0, 1 / 2), by norm_num [Q2]⟩⟩
instance : Nonempty Q1 := ⟨⟨0, by norm_num [Q1]⟩⟩

noncomputable def GchQ (b : ℝ → ℝ) : ℝ :=
  max 0 (⨆ q : Q3, b (q.1.2.1 : ℝ) - chord b (q.1.1 : ℝ) (q.1.2.2 : ℝ) (q.1.2.1 : ℝ))
noncomputable def LamQ (y : ℝ) (b : ℝ → ℝ) : ℝ :=
  max (b 1 - y) (⨆ q : Q2, endTerm b y (q.1.1 : ℝ) (q.1.2 : ℝ)) + y
noncomputable def XQ (κ : ℝ) (b : ℝ → ℝ) : ℝ := max 0 (⨆ q : Q1, b (1 - (q.1 : ℝ)) - κ * q.1)
noncomputable def oscQ (b : ℝ → ℝ) : ℝ := ⨆ q : Q1 × Q1, b (q.2.1 : ℝ) - b (q.1.1 : ℝ)

/-- `X(b) = max(0, sup_{d∈[0,1]} (b(1−d) − κd))` (DT11 §4 Lemma 3(b)). -/
noncomputable def Xr (κ : ℝ) (b : ℝ → ℝ) : ℝ :=
  max 0 (sSup {v | ∃ d ∈ Icc (0 : ℝ) 1, v = b (1 - d) - κ * d})
/-- Oscillation (range) `sup_{s,x∈[0,1]} (b x − b s)` (DT11 P3, Kuiper/Vervaat). -/
noncomputable def osc (b : ℝ → ℝ) : ℝ :=
  sSup {v | ∃ s ∈ Icc (0 : ℝ) 1, ∃ x ∈ Icc (0 : ℝ) 1, v = b x - b s}

lemma clipExt_cont {b : ℝ → ℝ} (hb : ContinuousOn b (Icc 0 1)) : Continuous (fun z => b (clip z)) :=
  hb.comp_continuous continuous_clip clip_mem

/-! ### Agreement on continuous paths -/

theorem GchQ_eq {b : ℝ → ℝ} (hb : ContinuousOn b (Icc 0 1)) : GchQ b = Gch b := by
  obtain ⟨C, hC⟩ := bound_of_cont hb
  have hbdd := gapSet_bdd hC
  have hsub : ∀ q : Q3, b (q.1.2.1 : ℝ) - chord b (q.1.1 : ℝ) (q.1.2.2 : ℝ) (q.1.2.1 : ℝ)
      ∈ gapSet b := fun q => ⟨_, _, _, q.2.1, q.2.2.1, q.2.2.2.1, q.2.2.2.2, rfl⟩
  have hbddR : BddAbove (range fun q : Q3 =>
      b (q.1.2.1 : ℝ) - chord b (q.1.1 : ℝ) (q.1.2.2 : ℝ) (q.1.2.1 : ℝ)) :=
    hbdd.mono (range_subset_iff.mpr hsub)
  unfold GchQ Gch
  congr 1
  apply le_antisymm
  · exact ciSup_le fun q => le_csSup hbdd (hsub q)
  · apply csSup_le (gapSet_nonempty b)
    rintro _ ⟨s, x, x', hs, hsx, hxx', hx', rfl⟩
    have he := clipExt_cont hb
    let G : ℝ × ℝ × ℝ → ℝ := fun p => b (clip p.2.1) -
      ((p.2.2 - p.2.1) / (p.2.2 - p.1) * b (clip p.1) + (p.2.1 - p.1) / (p.2.2 - p.1) * b (clip p.2.2))
    have hGI : ∀ p : ℝ × ℝ × ℝ, p.1 ∈ Icc (0 : ℝ) 1 → p.2.1 ∈ Icc (0 : ℝ) 1 →
        p.2.2 ∈ Icc (0 : ℝ) 1 → G p = b p.2.1 - chord b p.1 p.2.2 p.2.1 := by
      intro p h1 h2 h3
      simp only [G, chord, clip_of_mem h1, clip_of_mem h2, clip_of_mem h3]
    have hd : (x' - s) ≠ 0 := by linarith
    have hcont : ContinuousAt G (s, x, x') := by
      have c1 : Continuous fun p : ℝ × ℝ × ℝ => b (clip p.2.1) := he.comp (by fun_prop)
      have c0 : Continuous fun p : ℝ × ℝ × ℝ => b (clip p.1) := he.comp (by fun_prop)
      have c2 : Continuous fun p : ℝ × ℝ × ℝ => b (clip p.2.2) := he.comp (by fun_prop)
      have w1 : ContinuousAt (fun p : ℝ × ℝ × ℝ => (p.2.2 - p.2.1) / (p.2.2 - p.1)) (s, x, x') :=
        ContinuousAt.div (by fun_prop) (by fun_prop) hd
      have w2 : ContinuousAt (fun p : ℝ × ℝ × ℝ => (p.2.1 - p.1) / (p.2.2 - p.1)) (s, x, x') :=
        ContinuousAt.div (by fun_prop) (by fun_prop) hd
      exact c1.continuousAt.sub ((w1.mul c0.continuousAt).add (w2.mul c2.continuousAt))
    rw [← hGI (s, x, x') ⟨hs, by linarith⟩ ⟨by linarith, by linarith⟩ ⟨by linarith, hx'⟩]
    apply le_of_near hcont
    intro ρ hρ
    obtain ⟨q₁, q₂, q₃, h0, h12, h23, h3, d1, d2, d3⟩ := rat_triple hs hsx hxx' hx' hρ
    refine ⟨((q₁ : ℝ), (q₂ : ℝ), (q₃ : ℝ)), ?_, ?_⟩
    · simp only [Prod.dist_eq, Real.dist_eq]
      exact max_lt d1 (max_lt d2 d3)
    · rw [hGI _ ⟨h0, by linarith⟩ ⟨by linarith, by linarith⟩ ⟨by linarith, h3⟩]
      exact le_ciSup hbddR (⟨(q₁, q₂, q₃), h0, h12, h23, h3⟩ : Q3)

theorem LamQ_eq {b : ℝ → ℝ} (hb : ContinuousOn b (Icc 0 1)) {y : ℝ} (hy : 0 < y) :
    LamQ y b = Lam b y := by
  have hbdd := endSet_bdd hb hy
  have hsub : ∀ q : Q2, endTerm b y (q.1.1 : ℝ) (q.1.2 : ℝ) ∈ endSet b y :=
    fun q => ⟨_, _, q.2.1, q.2.2.1, q.2.2.2, rfl⟩
  have hbddR : BddAbove (range fun q : Q2 => endTerm b y (q.1.1 : ℝ) (q.1.2 : ℝ)) :=
    hbdd.mono (range_subset_iff.mpr hsub)
  unfold LamQ Lam Lhat
  congr 2
  apply le_antisymm
  · exact ciSup_le fun q => le_csSup hbdd (hsub q)
  · apply csSup_le (endSet_nonempty b y)
    rintro _ ⟨s, x, hs, hsx, hx1, rfl⟩
    have he := clipExt_cont hb
    let G : ℝ × ℝ → ℝ := fun p =>
      ((1 - p.1) * (b (clip p.2) - y) - (1 - p.2) * (b (clip p.1) + y)) / (p.2 - p.1)
    have hGI : ∀ p : ℝ × ℝ, p.1 ∈ Icc (0 : ℝ) 1 → p.2 ∈ Icc (0 : ℝ) 1 →
        G p = endTerm b y p.1 p.2 := by
      intro p h1 h2
      simp only [G, endTerm, clip_of_mem h1, clip_of_mem h2]
    have hd : (x - s) ≠ 0 := by linarith
    have hcont : ContinuousAt G (s, x) := by
      have c1 : Continuous fun p : ℝ × ℝ => b (clip p.2) := he.comp (by fun_prop)
      have c0 : Continuous fun p : ℝ × ℝ => b (clip p.1) := he.comp (by fun_prop)
      have num : Continuous fun p : ℝ × ℝ =>
          (1 - p.1) * (b (clip p.2) - y) - (1 - p.2) * (b (clip p.1) + y) :=
        ((continuous_const.sub continuous_fst).mul (c1.sub continuous_const)).sub
          ((continuous_const.sub continuous_snd).mul (c0.add continuous_const))
      exact ContinuousAt.div num.continuousAt (by fun_prop) hd
    rw [← hGI (s, x) ⟨hs, by linarith⟩ ⟨by linarith, hx1.le⟩]
    apply le_of_near hcont
    intro ρ hρ
    obtain ⟨q₁, q₂, h0, h12, h2, d1, d2⟩ := rat_pair hs hsx hx1 hρ
    refine ⟨((q₁ : ℝ), (q₂ : ℝ)), ?_, ?_⟩
    · simp only [Prod.dist_eq, Real.dist_eq]
      exact max_lt d1 d2
    · rw [hGI _ ⟨h0, by linarith⟩ ⟨by linarith, h2.le⟩]
      exact le_ciSup hbddR (⟨(q₁, q₂), h0, h12, h2⟩ : Q2)

theorem XQ_eq {b : ℝ → ℝ} (hb : ContinuousOn b (Icc 0 1)) (κ : ℝ) : XQ κ b = Xr κ b := by
  obtain ⟨C, hC⟩ := bound_of_cont hb
  have hbdd : BddAbove {v | ∃ d ∈ Icc (0 : ℝ) 1, v = b (1 - d) - κ * d} := by
    refine ⟨C + |κ|, ?_⟩
    rintro _ ⟨d, hd, rfl⟩
    have h1 := le_abs_self (b (1 - d))
    have h2 := hC (1 - d) ⟨by linarith [hd.2], by linarith [hd.1]⟩
    have h3 : -(κ * d) ≤ |κ| := by
      have := neg_abs_le κ
      have : |κ| * d ≤ |κ| := by nlinarith [abs_nonneg κ, hd.1, hd.2]
      nlinarith [hd.1, neg_abs_le κ, le_abs_self κ]
    linarith
  have hsub : ∀ q : Q1, b (1 - (q.1 : ℝ)) - κ * q.1 ∈ {v | ∃ d ∈ Icc (0 : ℝ) 1, v = b (1 - d) - κ * d} :=
    fun q => ⟨q.1, q.2, rfl⟩
  have hbddR : BddAbove (range fun q : Q1 => b (1 - (q.1 : ℝ)) - κ * q.1) :=
    hbdd.mono (range_subset_iff.mpr hsub)
  unfold XQ Xr
  congr 1
  apply le_antisymm
  · exact ciSup_le fun q => le_csSup hbdd (hsub q)
  · apply csSup_le (show ({v | ∃ d ∈ Icc (0 : ℝ) 1, v = b (1 - d) - κ * d} : Set ℝ).Nonempty
      from ⟨_, 0, ⟨le_rfl, zero_le_one⟩, rfl⟩)
    rintro _ ⟨d, hd, rfl⟩
    have he := clipExt_cont hb
    let G : ℝ → ℝ := fun p => b (clip (1 - p)) - κ * p
    have hGI : ∀ p ∈ Icc (0 : ℝ) 1, G p = b (1 - p) - κ * p := fun p hp => by
      simp only [G, clip_of_mem (show 1 - p ∈ Icc (0 : ℝ) 1 from ⟨by linarith [hp.2], by linarith [hp.1]⟩)]
    have hcont : ContinuousAt G d := (Continuous.sub (he.comp (by fun_prop)) (by fun_prop)).continuousAt
    rw [← hGI d hd]
    apply le_of_near hcont
    intro ρ hρ
    obtain ⟨q, hq, hqd⟩ := rat_single hd hρ
    refine ⟨q, by rw [Real.dist_eq]; exact hqd, ?_⟩
    rw [hGI _ hq]
    exact le_ciSup hbddR (⟨q, hq⟩ : Q1)

theorem oscQ_eq {b : ℝ → ℝ} (hb : ContinuousOn b (Icc 0 1)) : oscQ b = osc b := by
  obtain ⟨C, hC⟩ := bound_of_cont hb
  have hbdd : BddAbove {v | ∃ s ∈ Icc (0 : ℝ) 1, ∃ x ∈ Icc (0 : ℝ) 1, v = b x - b s} := by
    refine ⟨2 * C, ?_⟩
    rintro _ ⟨s, hs, x, hx, rfl⟩
    have := abs_le.mp (hC s hs); have := abs_le.mp (hC x hx); linarith
  have hsub : ∀ q : Q1 × Q1, b (q.2.1 : ℝ) - b (q.1.1 : ℝ) ∈
      {v | ∃ s ∈ Icc (0 : ℝ) 1, ∃ x ∈ Icc (0 : ℝ) 1, v = b x - b s} :=
    fun q => ⟨q.1.1, q.1.2, q.2.1, q.2.2, rfl⟩
  have hbddR : BddAbove (range fun q : Q1 × Q1 => b (q.2.1 : ℝ) - b (q.1.1 : ℝ)) :=
    hbdd.mono (range_subset_iff.mpr hsub)
  unfold oscQ osc
  apply le_antisymm
  · exact ciSup_le fun q => le_csSup hbdd (hsub q)
  · apply csSup_le (show ({v | ∃ s ∈ Icc (0 : ℝ) 1, ∃ x ∈ Icc (0 : ℝ) 1, v = b x - b s} : Set ℝ).Nonempty
      from ⟨_, 0, ⟨le_rfl, zero_le_one⟩, 0, ⟨le_rfl, zero_le_one⟩, rfl⟩)
    rintro _ ⟨s, hs, x, hx, rfl⟩
    have he := clipExt_cont hb
    let G : ℝ × ℝ → ℝ := fun p => b (clip p.2) - b (clip p.1)
    have hGI : ∀ p : ℝ × ℝ, p.1 ∈ Icc (0 : ℝ) 1 → p.2 ∈ Icc (0 : ℝ) 1 → G p = b p.2 - b p.1 :=
      fun p h1 h2 => by simp only [G, clip_of_mem h1, clip_of_mem h2]
    have hcont : ContinuousAt G (s, x) :=
      ((he.comp (by fun_prop)).sub (he.comp (by fun_prop))).continuousAt
    rw [← hGI (s, x) hs hx]
    apply le_of_near hcont
    intro ρ hρ
    obtain ⟨q₁, hq₁, d1⟩ := rat_single hs hρ
    obtain ⟨q₂, hq₂, d2⟩ := rat_single hx hρ
    refine ⟨((q₁ : ℝ), (q₂ : ℝ)), ?_, ?_⟩
    · simp only [Prod.dist_eq, Real.dist_eq]; exact max_lt d1 d2
    · rw [hGI _ hq₁ hq₂]
      exact le_ciSup hbddR ((⟨q₁, hq₁⟩, ⟨q₂, hq₂⟩) : Q1 × Q1)

/-- The oscillation bounds every increment. -/
lemma incr_le_osc {b : ℝ → ℝ} (hb : ContinuousOn b (Icc 0 1)) {s x : ℝ}
    (hs : s ∈ Icc (0 : ℝ) 1) (hx : x ∈ Icc (0 : ℝ) 1) : b x - b s ≤ osc b := by
  obtain ⟨C, hC⟩ := bound_of_cont hb
  refine le_csSup (s := {v | ∃ s ∈ Icc (0 : ℝ) 1, ∃ x ∈ Icc (0 : ℝ) 1, v = b x - b s})
    ⟨2 * C, ?_⟩ ⟨s, hs, x, hx, rfl⟩
  rintro _ ⟨s', hs', x', hx', rfl⟩
  have := abs_le.mp (hC s' hs'); have := abs_le.mp (hC x' hx'); linarith

lemma le_Xr {b : ℝ → ℝ} (hb : ContinuousOn b (Icc 0 1)) (κ : ℝ) {d : ℝ}
    (hd : d ∈ Icc (0 : ℝ) 1) : b (1 - d) - κ * d ≤ Xr κ b := by
  obtain ⟨C, hC⟩ := bound_of_cont hb
  refine (le_csSup (s := {v | ∃ d ∈ Icc (0 : ℝ) 1, v = b (1 - d) - κ * d})
    ⟨C + |κ|, ?_⟩ ⟨d, hd, rfl⟩).trans (le_max_right _ _)
  rintro _ ⟨d', hd', rfl⟩
  have h1 := le_abs_self (b (1 - d'))
  have h2 := hC (1 - d') ⟨by linarith [hd'.2], by linarith [hd'.1]⟩
  nlinarith [hd'.1, hd'.2, neg_abs_le κ, le_abs_self κ, abs_nonneg κ]

lemma Xr_nonneg (κ : ℝ) (b : ℝ → ℝ) : 0 ≤ Xr κ b := le_max_left _ _

lemma Xr_le {b : ℝ → ℝ} {κ w₀ : ℝ} (hκ : 0 ≤ κ) (hb1 : b 1 = 0)
    (hw : ∀ s ∈ Icc (0 : ℝ) 1, ∀ x ∈ Icc (0 : ℝ) 1, b x - b s ≤ w₀)
    {x : ℝ} (hx : x ∈ Icc (0 : ℝ) 1) : Xr κ b ≤ b x + w₀ := by
  have h1I : (1 : ℝ) ∈ Icc (0 : ℝ) 1 := ⟨zero_le_one, le_rfl⟩
  apply max_le
  · have := hw x hx 1 h1I; rw [hb1] at this; linarith
  · apply csSup_le (show ({v | ∃ d ∈ Icc (0 : ℝ) 1, v = b (1 - d) - κ * d} : Set ℝ).Nonempty
      from ⟨_, 0, ⟨le_rfl, zero_le_one⟩, rfl⟩)
    rintro _ ⟨d, hd, rfl⟩
    have := hw x hx (1 - d) ⟨by linarith [hd.2], by linarith [hd.1]⟩
    nlinarith [hd.1]

/-! ### Measurability on (ℝ → ℝ) with the product σ-algebra -/

theorem measurable_GchQ : Measurable GchQ := by
  have h : ∀ q : Q3, Measurable fun b : ℝ → ℝ =>
      b (q.1.2.1 : ℝ) - chord b (q.1.1 : ℝ) (q.1.2.2 : ℝ) (q.1.2.1 : ℝ) := fun q => by
    unfold chord; fun_prop
  exact measurable_const.max (Measurable.iSup h)

theorem measurable_LamQ (y : ℝ) : Measurable (LamQ y) := by
  have h : ∀ q : Q2, Measurable fun b : ℝ → ℝ => endTerm b y (q.1.1 : ℝ) (q.1.2 : ℝ) :=
    fun q => by unfold endTerm; fun_prop
  have h1 : Measurable fun b : ℝ → ℝ => b 1 - y := by fun_prop
  exact (h1.max (Measurable.iSup h)).add_const y

theorem measurable_XQ (κ : ℝ) : Measurable (XQ κ) := by
  have h : ∀ q : Q1, Measurable fun b : ℝ → ℝ => b (1 - (q.1 : ℝ)) - κ * q.1 := fun q => by
    fun_prop
  exact measurable_const.max (Measurable.iSup h)

theorem measurable_oscQ : Measurable oscQ := by
  have h : ∀ q : Q1 × Q1, Measurable fun b : ℝ → ℝ => b (q.2.1 : ℝ) - b (q.1.1 : ℝ) :=
    fun q => by fun_prop
  exact Measurable.iSup h

end Seam

namespace Seam

open MeasureTheory ProbabilityTheory

/-! ## Module C3: Lemma B as a probability statement (DT10 §6) -/

variable {Ω : Type*} [MeasurableSpace Ω]

lemma Lam_pieceL (B : ℝ → ℝ) {u : ℝ} (hu0 : 0 < u) (t : ℝ) :
    Lam (pieceL B u) t = √u * Lam (bL B u) (t / √u) := by
  have ha : 0 < √u := Real.sqrt_pos.mpr hu0
  rw [Lam_congr (pieceL_eq B hu0), ← Lam_smul _ ha]; congr 1; field_simp

lemma measurable_bL {B : Ω → ℝ → ℝ} (hBm : ∀ z, Measurable fun ω => B ω z) (u : ℝ) :
    Measurable fun ω => bL (B ω) u := by
  rw [measurable_pi_iff]; intro s
  unfold bL pieceL
  exact ((hBm _).sub ((hBm u).const_mul _)).div_const _

lemma measurable_bR {B : Ω → ℝ → ℝ} (hBm : ∀ z, Measurable fun ω => B ω z) (u : ℝ) :
    Measurable fun ω => bR (B ω) u := by
  rw [measurable_pi_iff]; intro s
  unfold bR pieceR
  exact (((hBm u).const_mul _).sub (hBm _)).div_const _

section meas
variable {B : Ω → ℝ → ℝ} (hBc : ∀ ω, ContinuousOn (B ω) (Icc 0 1))
  (hBm : ∀ z, Measurable fun ω => B ω z) {u : ℝ} (hu0 : 0 < u) (hu1 : u < 1)
include hBc hBm hu0 hu1

lemma measurable_Gch_bL : Measurable fun ω => Gch (bL (B ω) u) := by
  have : (fun ω => Gch (bL (B ω) u)) = GchQ ∘ fun ω => bL (B ω) u :=
    funext fun ω => (GchQ_eq (bL_cont (hBc ω) hu0 hu1)).symm
  rw [this]; exact measurable_GchQ.comp (measurable_bL hBm u)

lemma measurable_Gch_bR : Measurable fun ω => Gch (bR (B ω) u) := by
  have : (fun ω => Gch (bR (B ω) u)) = GchQ ∘ fun ω => bR (B ω) u :=
    funext fun ω => (GchQ_eq (bR_cont (hBc ω) hu0 hu1)).symm
  rw [this]; exact measurable_GchQ.comp (measurable_bR hBm u)

lemma measurable_Lam_bL {y : ℝ} (hy : 0 < y) : Measurable fun ω => Lam (bL (B ω) u) y := by
  have : (fun ω => Lam (bL (B ω) u) y) = LamQ y ∘ fun ω => bL (B ω) u :=
    funext fun ω => (LamQ_eq (bL_cont (hBc ω) hu0 hu1) hy).symm
  rw [this]; exact (measurable_LamQ y).comp (measurable_bL hBm u)

lemma measurable_Lam_bR {y : ℝ} (hy : 0 < y) : Measurable fun ω => Lam (bR (B ω) u) y := by
  have : (fun ω => Lam (bR (B ω) u) y) = LamQ y ∘ fun ω => bR (B ω) u :=
    funext fun ω => (LamQ_eq (bR_cont (hBc ω) hu0 hu1) hy).symm
  rw [this]; exact (measurable_LamQ y).comp (measurable_bR hBm u)

lemma measurable_Xr_bR (κ : ℝ) : Measurable fun ω => Xr κ (bR (B ω) u) := by
  have : (fun ω => Xr κ (bR (B ω) u)) = XQ κ ∘ fun ω => bR (B ω) u :=
    funext fun ω => (XQ_eq (bR_cont (hBc ω) hu0 hu1) κ).symm
  rw [this]; exact (measurable_XQ κ).comp (measurable_bR hBm u)

lemma measurable_osc_bR : Measurable fun ω => osc (bR (B ω) u) := by
  have : (fun ω => osc (bR (B ω) u)) = oscQ ∘ fun ω => bR (B ω) u :=
    funext fun ω => (oscQ_eq (bR_cont (hBc ω) hu0 hu1)).symm
  rw [this]; exact measurable_oscQ.comp (measurable_bR hBm u)

omit hBm [MeasurableSpace Ω] in
/-- The seam criterion in the normalised pieces (`r > 0`). -/
lemma seamFit_iff_b {r : ℝ} (hr : 0 < r) (ω : Ω) :
    SeamFit (B ω) u r ↔ √u * Gch (bL (B ω) u) ≤ 2 * r ∧ √(1 - u) * Gch (bR (B ω) u) ≤ 2 * r ∧
      (√u * Lam (bL (B ω) u) (r / √u) - r) + (√(1 - u) * Lam (bR (B ω) u) (r / √(1 - u)) - r)
        ≤ 0 := by
  rw [seamFit_iff (hBc ω) hu0 hu1 hr, Gch_pieceL _ hu0, Gch_pieceR _ hu1]
  have e1 : Lhat (pieceL (B ω) u) r = √u * Lam (bL (B ω) u) (r / √u) - r := by
    rw [← Lam_pieceL _ hu0]; unfold Lam; ring
  have e2 : Lhat (pieceR (B ω) u) r = √(1 - u) * Lam (bR (B ω) u) (r / √(1 - u)) - r := by
    rw [← Lam_pieceR _ hu1]; unfold Lam; ring
  rw [e1, e2]

/-- **The event `{d_u ≥ t}` is measurable**: it is `⋂_{r ∈ ℚ, 0<r<t} {no fit at r}`. -/
theorem measurableSet_dK_ge (t : ℝ) : MeasurableSet {ω | t ≤ dK u (B ω)} := by
  have key : ∀ ω, t ≤ dK u (B ω) ↔
      ∀ r ∈ {r : ℚ | 0 < (r : ℝ) ∧ (r : ℝ) < t}, ¬ SeamFit (B ω) u r := by
    intro ω
    constructor
    · rintro h r ⟨hr0, hrt⟩ hf
      have := csInf_le (dK_set_bdd (B ω) u) ⟨hr0.le, hf⟩
      unfold dK at h; linarith
    · intro h
      by_contra hlt
      push Not at hlt
      obtain ⟨r, hr1, hr2⟩ := exists_rat_btwn hlt
      have hr0 : 0 < (r : ℝ) := lt_of_le_of_lt (dK_nonneg _ _) hr1
      exact h r ⟨hr0, hr2⟩ (seamFit_of_dK_lt (hBc ω) hr1)
  have hset : {ω | t ≤ dK u (B ω)} =
      ⋂ r ∈ {r : ℚ | 0 < (r : ℝ) ∧ (r : ℝ) < t}, {ω | ¬ SeamFit (B ω) u r} := by
    ext ω; simp only [mem_ofPred_eq, mem_iInter]; exact key ω
  rw [hset]
  refine MeasurableSet.biInter (Set.to_countable _) (fun r hr => ?_)
  have hr0 : 0 < (r : ℝ) := hr.1
  have hm : MeasurableSet {ω | √u * Gch (bL (B ω) u) ≤ 2 * r ∧
      √(1 - u) * Gch (bR (B ω) u) ≤ 2 * r ∧
      (√u * Lam (bL (B ω) u) (r / √u) - r) + (√(1 - u) * Lam (bR (B ω) u) (r / √(1 - u)) - r)
        ≤ 0} := by
    have m1 := measurable_Gch_bL hBc hBm hu0 hu1
    have m2 := measurable_Gch_bR hBc hBm hu0 hu1
    have m3 := measurable_Lam_bL hBc hBm hu0 hu1 (y := (r : ℝ) / √u)
      (div_pos hr0 (Real.sqrt_pos.mpr hu0))
    have m4 := measurable_Lam_bR hBc hBm hu0 hu1 (y := (r : ℝ) / √(1 - u))
      (div_pos hr0 (Real.sqrt_pos.mpr (by linarith)))
    exact (measurableSet_le (m1.const_mul _) measurable_const).inter
      ((measurableSet_le (m2.const_mul _) measurable_const).inter
        (measurableSet_le (((m3.const_mul _).sub_const _).add ((m4.const_mul _).sub_const _))
          measurable_const))
  have : {ω | ¬ SeamFit (B ω) u r} = {ω | √u * Gch (bL (B ω) u) ≤ 2 * r ∧
      √(1 - u) * Gch (bR (B ω) u) ≤ 2 * r ∧
      (√u * Lam (bL (B ω) u) (r / √u) - r) + (√(1 - u) * Lam (bR (B ω) u) (r / √(1 - u)) - r)
        ≤ 0}ᶜ := by
    ext ω; simp only [mem_ofPred_eq, mem_compl_iff]; rw [seamFit_iff_b hBc hu0 hu1 hr0 ω]
  rw [this]; exact hm.compl

end meas

/-- Upper tail of a variable with a continuous CDF `F`: `μ{z ≤ G} ≤ 1 − F z`. -/
lemma tail_le {μ : Measure Ω} [IsProbabilityMeasure μ] {G : Ω → ℝ} (hGm : Measurable G)
    {F : ℝ → ℝ} (hF01 : ∀ z, 0 ≤ F z ∧ F z ≤ 1) (hFc : Continuous F)
    (hlaw : ∀ z, μ {ω | G ω ≤ z} = ENNReal.ofReal (F z)) (z : ℝ) :
    μ {ω | z ≤ G ω} ≤ ENNReal.ofReal (1 - F z) := by
  have hlt : MeasurableSet {ω | G ω < z} := measurableSet_lt hGm measurable_const
  have hc : {ω | z ≤ G ω} = {ω | G ω < z}ᶜ := by ext ω; simp
  rw [hc, prob_compl_eq_one_sub hlt]
  have hlow : ENNReal.ofReal (F z) ≤ μ {ω | G ω < z} := by
    have ht : Filter.Tendsto (fun n : ℕ => z - 1 / ((n : ℝ) + 1)) Filter.atTop (nhds z) := by
      have := (tendsto_one_div_add_atTop_nhds_zero_nat (𝕜 := ℝ)).const_sub z
      simpa using this
    have hT : Filter.Tendsto (fun n : ℕ => ENNReal.ofReal (F (z - 1 / ((n : ℝ) + 1))))
        Filter.atTop (nhds (ENNReal.ofReal (F z))) :=
      (ENNReal.continuous_ofReal.tendsto _).comp ((hFc.tendsto z).comp ht)
    refine le_of_tendsto' hT (fun n => ?_)
    rw [← hlaw]
    apply measure_mono
    intro ω hω
    simp only [mem_ofPred_eq] at hω ⊢
    have : 0 < 1 / ((n : ℝ) + 1) := by positivity
    linarith
  calc 1 - μ {ω | G ω < z} ≤ 1 - ENNReal.ofReal (F z) := tsub_le_tsub_left hlow 1
    _ = ENNReal.ofReal (1 - F z) := by
      rw [ENNReal.ofReal_sub 1 (hF01 z).1, ENNReal.ofReal_one]

/-- Lower bound of an upper tail from the CDF: `1 − F z ≤ μ{z ≤ G}` (no continuity needed). -/
lemma tail_ge {μ : Measure Ω} [IsProbabilityMeasure μ] {G : Ω → ℝ} (hGm : Measurable G)
    {F : ℝ → ℝ} (hF01 : ∀ z, 0 ≤ F z ∧ F z ≤ 1)
    (hlaw : ∀ z, μ {ω | G ω ≤ z} = ENNReal.ofReal (F z)) (z : ℝ) :
    1 - F z ≤ (μ {ω | z ≤ G ω}).toReal := by
  have hle : MeasurableSet {ω | G ω ≤ z} := measurableSet_le hGm measurable_const
  have h1 : μ {ω | z < G ω} = ENNReal.ofReal (1 - F z) := by
    have hc : {ω | z < G ω} = {ω | G ω ≤ z}ᶜ := by ext ω; simp
    rw [hc, prob_compl_eq_one_sub hle, hlaw, ENNReal.ofReal_sub 1 (hF01 z).1, ENNReal.ofReal_one]
  have h2 : μ {ω | z < G ω} ≤ μ {ω | z ≤ G ω} :=
    measure_mono (fun ω (h : z < G ω) => (show z ≤ G ω from h.le))
  have h3 := ENNReal.toReal_mono (measure_ne_top μ _) h2
  rw [h1, ENNReal.toReal_ofReal (by linarith [(hF01 z).2])] at h3
  exact h3

/-- **Lemma B** (DT10 §6) as a probability inequality, for the non-strict tail `P(d_u ≥ t)`.
Trusted inputs: (S) `b₁ ⟂ b₂` (Lemma S(i), Gaussian splitting of the bridge), as maps into
`ℝ → ℝ` (clipped pieces); (M) the CDF of `G(b₁)` is the continuous `F` (BP11 + sign symmetry). -/
theorem lemmaB_prob (μ : Measure Ω) [IsProbabilityMeasure μ] {B : Ω → ℝ → ℝ}
    (hBc : ∀ ω, ContinuousOn (B ω) (Icc 0 1)) (hBm : ∀ z, Measurable fun ω => B ω z)
    {u t : ℝ} (hu0 : 0 < u) (hu1 : u < 1) (ht : 0 < t)
    {F : ℝ → ℝ} (hF01 : ∀ z, 0 ≤ F z ∧ F z ≤ 1) (hFc : Continuous F)
    (hind : IndepFun (fun ω => bL (B ω) u) (fun ω => bR (B ω) u) μ)
    (hlaw : ∀ z, μ {ω | Gch (bL (B ω) u) ≤ z} = ENNReal.ofReal (F z)) :
    μ {ω | t ≤ dK u (B ω)} ≤ μ {ω | 2 * t ≤ √(1 - u) * Gch (bR (B ω) u)} +
      ∫⁻ ω, (if √(1 - u) * Gch (bR (B ω) u) < 2 * t then
        ENNReal.ofReal (1 - F ((2 * t - √(1 - u) * Lam (bR (B ω) u) (t / √(1 - u))) / √u))
        else 0) ∂μ := by
  set ε := √(1 - u) with hε
  set a := √u with ha
  have hε0 : 0 < ε := Real.sqrt_pos.mpr (by linarith)
  have ha0 : 0 < a := Real.sqrt_pos.mpr hu0
  set y := t / ε with hy
  have hy0 : 0 < y := div_pos ht hε0
  set G1 : Ω → ℝ := fun ω => Gch (bL (B ω) u) with hG1
  set Y : Ω → ℝ × ℝ := fun ω => (Gch (bR (B ω) u), Lam (bR (B ω) u) y) with hY
  have hG1m : Measurable G1 := measurable_Gch_bL hBc hBm hu0 hu1
  have hYm : Measurable Y :=
    (measurable_Gch_bR hBc hBm hu0 hu1).prodMk (measurable_Lam_bR hBc hBm hu0 hu1 hy0)
  -- independence of the two functionals (push-forward of Lemma S(i))
  have hindGY : IndepFun G1 Y μ := by
    have e1 : G1 = GchQ ∘ fun ω => bL (B ω) u :=
      funext fun ω => (GchQ_eq (bL_cont (hBc ω) hu0 hu1)).symm
    have e2 : Y = (fun f => (GchQ f, LamQ y f)) ∘ fun ω => bR (B ω) u :=
      funext fun ω => by
        simp only [hY, Function.comp, GchQ_eq (bR_cont (hBc ω) hu0 hu1),
          LamQ_eq (bR_cont (hBc ω) hu0 hu1) hy0]
    rw [e1, e2]
    exact hind.comp measurable_GchQ (measurable_GchQ.prodMk (measurable_LamQ y))
  -- deterministic inclusion (Lemma B core)
  set E' : Set Ω := {ω | ε * (Y ω).1 < 2 * t ∧ (2 * t - ε * (Y ω).2) / a ≤ G1 ω} with hE'
  have hincl : {ω | t ≤ dK u (B ω)} ⊆ {ω | 2 * t ≤ ε * Gch (bR (B ω) u)} ∪ E' := by
    intro ω hω
    rcases lemmaB_det (hBc ω) hu0 hu1 ht hω with h | h
    · exact Or.inl h
    · exact Or.inr h
  -- the product-measure computation for E'
  set S : Set (ℝ × (ℝ × ℝ)) := {p | ε * p.2.1 < 2 * t ∧ (2 * t - ε * p.2.2) / a ≤ p.1} with hS
  have hSm : MeasurableSet S := by
    have h1 : Measurable fun p : ℝ × (ℝ × ℝ) => ε * p.2.1 := by fun_prop
    have h2 : Measurable fun p : ℝ × (ℝ × ℝ) => (2 * t - ε * p.2.2) / a := by fun_prop
    exact (measurableSet_lt h1 measurable_const).inter (measurableSet_le h2 measurable_fst)
  set g : ℝ × ℝ → ENNReal := fun q => if ε * q.1 < 2 * t then
      ENNReal.ofReal (1 - F ((2 * t - ε * q.2) / a)) else 0 with hg
  have hgm : Measurable g := by
    refine Measurable.ite (measurableSet_lt (by fun_prop) measurable_const) ?_ measurable_const
    exact ENNReal.measurable_ofReal.comp (measurable_const.sub (hFc.measurable.comp (by fun_prop)))
  have hE'eq : E' = (fun ω => (G1 ω, Y ω)) ⁻¹' S := rfl
  have hE'le : μ E' ≤ ∫⁻ ω, g (Y ω) ∂μ := by
    rw [hE'eq, ← Measure.map_apply (hG1m.prodMk hYm) hSm,
      (indepFun_iff_map_prod_eq_prod_map_map hG1m.aemeasurable hYm.aemeasurable).mp hindGY,
      Measure.prod_apply_symm hSm, ← lintegral_map hgm hYm]
    apply lintegral_mono
    intro q
    show (Measure.map G1 μ) ((fun x => (x, q)) ⁻¹' S) ≤ g q
    by_cases hq : ε * q.1 < 2 * t
    · have hpre : (fun x => (x, q)) ⁻¹' S = Ici ((2 * t - ε * q.2) / a) := by
        ext x; simp [hS, hq]
      rw [hpre, Measure.map_apply hG1m measurableSet_Ici]
      simp only [hg, hq, ↓reduceIte]
      exact tail_le hG1m hF01 hFc hlaw _
    · have hpre : (fun x => (x, q)) ⁻¹' S = ∅ := by
        ext x; simp [hS, hq]
      rw [hpre, measure_empty]
      exact zero_le
  calc μ {ω | t ≤ dK u (B ω)}
      ≤ μ ({ω | 2 * t ≤ ε * Gch (bR (B ω) u)} ∪ E') := measure_mono hincl
    _ ≤ μ {ω | 2 * t ≤ ε * Gch (bR (B ω) u)} + μ E' := measure_union_le _ _
    _ ≤ μ {ω | 2 * t ≤ ε * Gch (bR (B ω) u)} + ∫⁻ ω, g (Y ω) ∂μ := add_le_add le_rfl hE'le

end Seam

namespace Seam

open MeasureTheory ProbabilityTheory

/-! ## Module D: Lemma B joined to Theorem NE (`NearEdge.near_edge_band`) -/

variable {Ω : Type*} [MeasurableSpace Ω]

/-- Edge side: `1 − F(2t) ≤ P(d₁(B) ≥ t)` (uses `d₁ = G/2`, `dK_one`). -/
lemma edge_ge (μ : Measure Ω) [IsProbabilityMeasure μ] {B : Ω → ℝ → ℝ}
    (hBc : ∀ ω, ContinuousOn (B ω) (Icc 0 1)) (hBm : ∀ z, Measurable fun ω => B ω z)
    {F : ℝ → ℝ} (hF01 : ∀ z, 0 ≤ F z ∧ F z ≤ 1)
    (hlawB : ∀ z, μ {ω | Gch (B ω) ≤ z} = ENNReal.ofReal (F z)) (t : ℝ) :
    1 - F (2 * t) ≤ (μ {ω | t ≤ dK 1 (B ω)}).toReal := by
  have hGm : Measurable fun ω => Gch (B ω) := by
    have : (fun ω => Gch (B ω)) = GchQ ∘ B := funext fun ω => (GchQ_eq (hBc ω)).symm
    rw [this]; exact measurable_GchQ.comp (measurable_pi_iff.mpr hBm)
  have hset : {ω | t ≤ dK 1 (B ω)} = {ω | 2 * t ≤ Gch (B ω)} := by
    ext ω; simp only [mem_ofPred_eq, dK_one (hBc ω)]; constructor <;> intro h <;> linarith
  rw [hset]; exact tail_ge hGm hF01 hlawB (2 * t)

/-- Lemma B in the real-valued form used by `NearEdge.near_edge_band` (hypothesis `hB`). -/
theorem lemmaB_real (μ : Measure Ω) [IsProbabilityMeasure μ] {B : Ω → ℝ → ℝ}
    (hBc : ∀ ω, ContinuousOn (B ω) (Icc 0 1)) (hBm : ∀ z, Measurable fun ω => B ω z)
    {u t : ℝ} (hu0 : 0 < u) (hu1 : u < 1) (ht : 0 < t)
    {F : ℝ → ℝ} (hF01 : ∀ z, 0 ≤ F z ∧ F z ≤ 1) (hFc : Continuous F)
    (hind : IndepFun (fun ω => bL (B ω) u) (fun ω => bR (B ω) u) μ)
    (hlaw : ∀ z, μ {ω | Gch (bL (B ω) u) ≤ z} = ENNReal.ofReal (F z)) :
    (μ {ω | t ≤ dK u (B ω)}).toReal ≤ (μ {ω | 2 * t ≤ √(1 - u) * Gch (bR (B ω) u)}).toReal +
      ∫ ω, (1 - F (2 * t / √u - √(1 - u) / √u * Lam (bR (B ω) u) (t / √(1 - u)))) ∂μ := by
  have h := lemmaB_prob μ hBc hBm hu0 hu1 ht hF01 hFc hind hlaw
  have hε0 : 0 < √(1 - u) := Real.sqrt_pos.mpr (by linarith)
  have ha0 : 0 < √u := Real.sqrt_pos.mpr hu0
  set f : Ω → ℝ := fun ω => if √(1 - u) * Gch (bR (B ω) u) < 2 * t then
      1 - F ((2 * t - √(1 - u) * Lam (bR (B ω) u) (t / √(1 - u))) / √u) else 0 with hf
  set g : Ω → ℝ := fun ω =>
      1 - F (2 * t / √u - √(1 - u) / √u * Lam (bR (B ω) u) (t / √(1 - u))) with hg
  have hLm := measurable_Lam_bR hBc hBm hu0 hu1 (div_pos ht hε0)
  have hGm := measurable_Gch_bR hBc hBm hu0 hu1
  have hfm : Measurable f := by
    refine Measurable.ite (measurableSet_lt (hGm.const_mul _) measurable_const) ?_ measurable_const
    exact measurable_const.sub (hFc.measurable.comp (((measurable_const.sub
      (hLm.const_mul _)).div_const _)))
  have hgm : Measurable g :=
    measurable_const.sub (hFc.measurable.comp (measurable_const.sub (hLm.const_mul _)))
  have hf0 : ∀ ω, 0 ≤ f ω := fun ω => by
    simp only [hf]; split_ifs
    · linarith [(hF01 ((2 * t - √(1 - u) * Lam (bR (B ω) u) (t / √(1 - u))) / √u)).2]
    · exact le_rfl
  have hf1 : ∀ ω, f ω ≤ 1 := fun ω => by
    simp only [hf]; split_ifs
    · linarith [(hF01 ((2 * t - √(1 - u) * Lam (bR (B ω) u) (t / √(1 - u))) / √u)).1]
    · norm_num
  have hfg : ∀ ω, f ω ≤ g ω := fun ω => by
    simp only [hf, hg]
    have e : (2 * t - √(1 - u) * Lam (bR (B ω) u) (t / √(1 - u))) / √u
        = 2 * t / √u - √(1 - u) / √u * Lam (bR (B ω) u) (t / √(1 - u)) := by
      field_simp
    split_ifs
    · rw [e]
    · linarith [(hF01 (2 * t / √u - √(1 - u) / √u * Lam (bR (B ω) u) (t / √(1 - u)))).2]
  have hfi : Integrable f μ :=
    Integrable.of_mem_Icc 0 1 hfm.aemeasurable (ae_of_all _ fun ω => ⟨hf0 ω, hf1 ω⟩)
  have hgi : Integrable g μ :=
    Integrable.of_mem_Icc 0 1 hgm.aemeasurable (ae_of_all _ fun ω =>
      ⟨by simp only [hg]; linarith [(hF01 (2 * t / √u - √(1 - u) / √u *
          Lam (bR (B ω) u) (t / √(1 - u)))).2],
       by simp only [hg]; linarith [(hF01 (2 * t / √u - √(1 - u) / √u *
          Lam (bR (B ω) u) (t / √(1 - u)))).1]⟩)
  have hlin : (∫⁻ ω, (if √(1 - u) * Gch (bR (B ω) u) < 2 * t then
        ENNReal.ofReal (1 - F ((2 * t - √(1 - u) * Lam (bR (B ω) u) (t / √(1 - u))) / √u))
        else 0) ∂μ) = ∫⁻ ω, ENNReal.ofReal (f ω) ∂μ := by
    congr 1; funext ω; simp only [hf]; split_ifs <;> simp
  rw [hlin] at h
  have hI : (∫⁻ ω, ENNReal.ofReal (f ω) ∂μ).toReal = ∫ ω, f ω ∂μ :=
    (integral_eq_lintegral_of_nonneg_ae (ae_of_all _ hf0) hfm.aestronglyMeasurable).symm
  have hIfin : ∫⁻ ω, ENNReal.ofReal (f ω) ∂μ ≠ ⊤ := by
    have : ∫⁻ ω, ENNReal.ofReal (f ω) ∂μ ≤ ∫⁻ _, 1 ∂μ :=
      lintegral_mono fun ω => by rw [← ENNReal.ofReal_one]; exact ENNReal.ofReal_le_ofReal (hf1 ω)
    rw [lintegral_const, measure_univ, mul_one] at this
    exact ne_top_of_le_ne_top ENNReal.one_ne_top this
  have h2 := ENNReal.toReal_mono (ENNReal.add_ne_top.mpr ⟨measure_ne_top μ _, hIfin⟩) h
  rw [ENNReal.toReal_add (measure_ne_top μ _) hIfin, hI] at h2
  have h3 : ∫ ω, f ω ∂μ ≤ ∫ ω, g ω ∂μ := integral_mono hfi hgi hfg
  linarith

/-- **Near-edge band with Lemma B in Lean** (Theorem NE + DT10 Lemma B + Lemma F(d)).
For `t ∈ [0.72, 0.88]` and `ε ∈ (0, 0.15]`: `P(d_{1−ε²}(B) ≥ t) < P(d₁(B) ≥ t)`.
Lean-checked deduction from the trusted inputs:
* (S) Lemma S(i): the clipped normalised pieces `b₁ = bL`, `b₂ = bR` at `u = 1 − ε²` are
  independent (Gaussian splitting of the Brownian bridge);
* (M) BP11 + sign symmetry: `G(b₁)` and `G(B)` have the continuous CDF `F = F_M`;
* UP/LO for `F` (Lemma 1, `NearEdge.lemma1_increments`, from P1, P2, P5);
* (P3) Kuiper/Vervaat for the oscillation of `b₂`; (P4) linear-boundary crossing for `X(b₂)`;
* the paths are continuous on `[0,1]` and their evaluations are measurable. -/
theorem near_edge_seam (μ : Measure Ω) [IsProbabilityMeasure μ] (B : Ω → ℝ → ℝ)
    (hBc : ∀ ω, ContinuousOn (B ω) (Icc 0 1)) (hBm : ∀ z, Measurable fun ω => B ω z)
    (F : ℝ → ℝ) (hF01 : ∀ z, 0 ≤ F z ∧ F z ≤ 1) (hFmono : Monotone F) (hFc : Continuous F)
    (hUP : ∀ c d : ℝ, 1 ≤ c → c ≤ d → d ≤ 1.8 → F d - F c ≤ ∫ z in c..d, NearEdge.I1 z)
    (hLO : ∀ c d : ℝ, 1.44 ≤ c → c ≤ d → d ≤ 1.8 →
      (1 - 0.065) * ∫ z in c..d, NearEdge.I1 z ≤ F d - F c)
    (t ε : ℝ) (ht : 0.72 ≤ t ∧ t ≤ 0.88) (hε : 0 < ε ∧ ε ≤ 0.15)
    (hind : IndepFun (fun ω => bL (B ω) (1 - ε ^ 2)) (fun ω => bR (B ω) (1 - ε ^ 2)) μ)
    (hlaw1 : ∀ z, μ {ω | Gch (bL (B ω) (1 - ε ^ 2)) ≤ z} = ENNReal.ofReal (F z))
    (hlawB : ∀ z, μ {ω | Gch (B ω) ≤ z} = ENNReal.ofReal (F z))
    (hKuiper : (μ {ω | 0.6 * t / ε < osc (bR (B ω) (1 - ε ^ 2))}).toReal
      ≤ NearEdge.S (0.6 * t / ε))
    (hcross : ∀ r : ℝ, 0 < r → μ {ω | r < Xr (1.4 * t / ε) (bR (B ω) (1 - ε ^ 2))}
      ≤ ENNReal.ofReal (Real.exp (-2 * r * (r + 1.4 * t / ε)))) :
    (μ {ω | t ≤ dK (1 - ε ^ 2) (B ω)}).toReal < (μ {ω | t ≤ dK 1 (B ω)}).toReal := by
  obtain ⟨ht1, ht2⟩ := ht
  obtain ⟨hε0, hε1⟩ := hε
  have ht0 : 0 < t := by linarith
  have hε2 : ε ^ 2 ≤ 0.0225 := by nlinarith
  have hu0 : 0 < 1 - ε ^ 2 := by linarith
  have hu1 : 1 - ε ^ 2 < 1 := by have := pow_pos hε0 2; linarith
  have hsq : √(1 - (1 - ε ^ 2)) = ε := by
    rw [show 1 - (1 - ε ^ 2) = ε ^ 2 by ring]; exact Real.sqrt_sq hε0.le
  set κ := 1.4 * t / ε with hκ
  set w₀ := 0.6 * t / ε with hw₀
  have hκ0 : 0 ≤ κ := by rw [hκ]; positivity
  have hwk : w₀ + κ ≤ 2 * (t / ε) := by rw [hw₀, hκ]; field_simp; ring_nf; exact le_rfl
  set T : Ω → ℝ := fun ω => Lam (bR (B ω) (1 - ε ^ 2)) (t / ε) with hT
  set X : Ω → ℝ := fun ω => Xr κ (bR (B ω) (1 - ε ^ 2)) with hX
  set A : Set Ω := {ω | osc (bR (B ω) (1 - ε ^ 2)) ≤ w₀} with hA
  set E : Set Ω := {ω | 2 * t ≤ ε * Gch (bR (B ω) (1 - ε ^ 2))} with hE
  have hTm : Measurable T := measurable_Lam_bR hBc hBm hu0 hu1 (div_pos ht0 hε0)
  have hXm : Measurable X := measurable_Xr_bR hBc hBm hu0 hu1 κ
  have hAm : MeasurableSet A := measurableSet_le (measurable_osc_bR hBc hBm hu0 hu1) measurable_const
  have hw : ∀ ω ∈ A, ∀ s ∈ Icc (0 : ℝ) 1, ∀ x ∈ Icc (0 : ℝ) 1,
      bR (B ω) (1 - ε ^ 2) x - bR (B ω) (1 - ε ^ 2) s ≤ w₀ := fun ω hω s hs x hx =>
    (incr_le_osc (bR_cont (hBc ω) hu0 hu1) hs hx).trans hω
  have hTX : ∀ ω ∈ A, T ω ≤ X ω := by
    intro ω hω
    have hb := bR_cont (hBc ω) hu0 hu1
    refine Lam_le_X hκ0 hwk (fun x hx => ?_) (fun x hx => Xr_le hκ0 (bR_one _ _) (hw ω hω) hx)
    have := le_Xr hb κ (d := 1 - x) ⟨by linarith [hx.2], by linarith [hx.1]⟩
    rwa [show 1 - (1 - x) = x by ring] at this
  have hAc : Aᶜ = {ω | w₀ < osc (bR (B ω) (1 - ε ^ 2))} := by ext ω; simp [hA]
  have hAS : (μ Aᶜ).toReal ≤ NearEdge.S w₀ := by rw [hAc]; exact hKuiper
  have hEA : E ⊆ Aᶜ := by
    intro ω hωE hωA
    have hG := Gch_le_range (hw ω hωA)
    have : ε * Gch (bR (B ω) (1 - ε ^ 2)) ≤ ε * w₀ := mul_le_mul_of_nonneg_left hG hε0.le
    have e : ε * w₀ = 0.6 * t := by rw [hw₀]; field_simp
    simp only [hE, mem_ofPred_eq] at hωE
    linarith
  have hES : (μ E).toReal ≤ NearEdge.S w₀ :=
    (ENNReal.toReal_mono (measure_ne_top μ _) (measure_mono hEA)).trans hAS
  have hB := lemmaB_real μ hBc hBm hu0 hu1 ht0 hF01 hFc hind hlaw1
  rw [hsq] at hB
  have hband := NearEdge.near_edge_band μ F t ε _ ⟨ht1, ht2⟩ ⟨hε0, hε1⟩ hF01 hFmono hUP hLO
    T X A E hTm hXm hAm (fun ω => Xr_nonneg _ _) hTX hcross hAS hES hB
  exact lt_of_lt_of_le hband (edge_ge μ hBc hBm hF01 hlawB t)

/-- The same with Lemma 1 discharged from the published inputs P1, P2, P5
(`NearEdge.lemma1_increments`). -/
theorem near_edge_seam_published
    {Ω' : Type*} [MeasurableSpace Ω'] (ν : Measure Ω') [IsProbabilityMeasure ν]
    (L : ℕ → Ω' → ℝ) (hLm : ∀ j, Measurable (L j)) (hL : ∀ j ω, 0 < L j ω ∧ L j ω ≤ 1)
    (P pE : ℝ → ℝ) (hPmono : Monotone P) (hP01 : ∀ v, 0 ≤ P v ∧ P v ≤ 1)
    (hPS : ∀ v, 1 ≤ v → P v = 1 - NearEdge.S v)
    (hPd : ∀ v, 1 ≤ v → HasDerivAt P (pE v) v)
    (hpEc : ContinuousOn pE (Ici 1)) (hpEm : Measurable pE)
    (hP5a : ∀ h : ℝ → ENNReal, Measurable h →
      ∫⁻ ω, ∑' j, h (L j ω) ∂ν = ∫⁻ ℓ in Ioo (0 : ℝ) 1, ENNReal.ofReal (1 / ℓ) * h ℓ)
    (hP5b : ∀ h k : ℝ → ENNReal, Measurable h → Measurable k →
      ∫⁻ ω, ∑' j, ∑' i, (if i = j then 0 else h (L j ω) * k (L i ω)) ∂ν
        = ∫⁻ p in NearEdge.simplex2, ENNReal.ofReal (1 / (p.1 * p.2)) * (h p.1 * k p.2))
    (F : ℝ → ℝ) (hF01 : ∀ z, 0 ≤ F z ∧ F z ≤ 1) (hFmono : Monotone F) (hFc : Continuous F)
    (hP1 : ∀ z ∈ Icc (1 : ℝ) 1.8, Filter.Tendsto
      (fun n => ∫ ω, ∏ i ∈ Finset.range n, P (z / √(L i ω)) ∂ν) Filter.atTop (nhds (F z)))
    (μ : Measure Ω) [IsProbabilityMeasure μ] (B : Ω → ℝ → ℝ)
    (hBc : ∀ ω, ContinuousOn (B ω) (Icc 0 1)) (hBm : ∀ z, Measurable fun ω => B ω z)
    (t ε : ℝ) (ht : 0.72 ≤ t ∧ t ≤ 0.88) (hε : 0 < ε ∧ ε ≤ 0.15)
    (hind : IndepFun (fun ω => bL (B ω) (1 - ε ^ 2)) (fun ω => bR (B ω) (1 - ε ^ 2)) μ)
    (hlaw1 : ∀ z, μ {ω | Gch (bL (B ω) (1 - ε ^ 2)) ≤ z} = ENNReal.ofReal (F z))
    (hlawB : ∀ z, μ {ω | Gch (B ω) ≤ z} = ENNReal.ofReal (F z))
    (hKuiper : (μ {ω | 0.6 * t / ε < osc (bR (B ω) (1 - ε ^ 2))}).toReal
      ≤ NearEdge.S (0.6 * t / ε))
    (hcross : ∀ r : ℝ, 0 < r → μ {ω | r < Xr (1.4 * t / ε) (bR (B ω) (1 - ε ^ 2))}
      ≤ ENNReal.ofReal (Real.exp (-2 * r * (r + 1.4 * t / ε)))) :
    (μ {ω | t ≤ dK (1 - ε ^ 2) (B ω)}).toReal < (μ {ω | t ≤ dK 1 (B ω)}).toReal := by
  obtain ⟨hUP, hLO⟩ := NearEdge.lemma1_increments ν L hLm hL P pE hPmono hP01 hPS hPd hpEc hpEm
    hP5a hP5b F hP1
  exact near_edge_seam μ B hBc hBm F hF01 hFmono hFc hUP hLO t ε ht hε hind hlaw1 hlawB hKuiper
    hcross

end Seam
