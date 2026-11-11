/-
DT8 Step 10 (gap G3): machine checks for the deterministic core of DT7-UNIF Theorem U1, Step A
(domination lemma). Statements mirror DT7-UNIF.md:
* comp_convexOn_of_monotone : g convex nondecreasing on t, f convex on s, f maps s into t => g ∘ f convex.
* comp_concaveOn_of_monotone: the concave version.
* flatten_left_*  : the monotone repair φ̃ = φ ∘ (max · r) of a convex branch with minimiser r is convex,
                    nondecreasing, and keeps the error band against any nondecreasing Γ.
* flatten_right_* : the same for a concave branch with maximiser r', φ̃ = φ ∘ (min · r').
* transfer_error  : the error of φ̃ ∘ F against Γ ∘ F is at most the error of φ̃ against Γ.
H&H Theorem 1 (dip = distance to U*) is an external published theorem (EXTERNAL-THEOREMS E1) and is not
formalised here.
-/
import Mathlib

open Set

/-- Convex nondecreasing outer function of a convex inner function is convex. -/
theorem comp_convexOn_of_monotone {s t : Set ℝ} {f g : ℝ → ℝ}
    (hf : ConvexOn ℝ s f) (hg : ConvexOn ℝ t g) (hgm : MonotoneOn g t) (hmap : MapsTo f s t) :
    ConvexOn ℝ s (g ∘ f) := by
  refine ⟨hf.1, ?_⟩
  intro x hx y hy a b ha hb hab
  have hxy : a • x + b • y ∈ s := hf.1 hx hy ha hb hab
  have h1 : f (a • x + b • y) ≤ a • f x + b • f y := hf.2 hx hy ha hb hab
  have hfx : f x ∈ t := hmap hx
  have hfy : f y ∈ t := hmap hy
  have hmid : a • f x + b • f y ∈ t := hg.1 hfx hfy ha hb hab
  calc (g ∘ f) (a • x + b • y) = g (f (a • x + b • y)) := rfl
    _ ≤ g (a • f x + b • f y) := hgm (hmap hxy) hmid h1
    _ ≤ a • g (f x) + b • g (f y) := hg.2 hfx hfy ha hb hab

/-- Concave nondecreasing outer function of a concave inner function is concave. -/
theorem comp_concaveOn_of_monotone {s t : Set ℝ} {f g : ℝ → ℝ}
    (hf : ConcaveOn ℝ s f) (hg : ConcaveOn ℝ t g) (hgm : MonotoneOn g t) (hmap : MapsTo f s t) :
    ConcaveOn ℝ s (g ∘ f) := by
  refine ⟨hf.1, ?_⟩
  intro x hx y hy a b ha hb hab
  have hxy : a • x + b • y ∈ s := hf.1 hx hy ha hb hab
  have h1 : a • f x + b • f y ≤ f (a • x + b • y) := hf.2 hx hy ha hb hab
  have hfx : f x ∈ t := hmap hx
  have hfy : f y ∈ t := hmap hy
  have hmid : a • f x + b • f y ∈ t := hg.1 hfx hfy ha hb hab
  calc a • g (f x) + b • g (f y) ≤ g (a • f x + b • f y) := hg.2 hfx hfy ha hb hab
    _ ≤ g (f (a • x + b • y)) := hgm hmid (hmap hxy) h1
    _ = (g ∘ f) (a • x + b • y) := rfl

/-- A convex function on [lo, hi] with a minimum at r is nondecreasing on [r, hi]. -/
theorem convex_monotone_after_min {φ : ℝ → ℝ} {lo hi r : ℝ}
    (hφ : ConvexOn ℝ (Icc lo hi) φ) (hr : r ∈ Icc lo hi) (hmin : ∀ z ∈ Icc lo hi, φ r ≤ φ z) :
    MonotoneOn φ (Icc r hi) := by
  intro x hx y hy hxy
  have hxs : x ∈ Icc lo hi := ⟨le_trans hr.1 hx.1, hx.2⟩
  have hys : y ∈ Icc lo hi := ⟨le_trans hr.1 hy.1, hy.2⟩
  have hseg : x ∈ segment ℝ r y := by
    rw [segment_eq_Icc (le_trans hx.1 hxy)]
    exact ⟨hx.1, hxy⟩
  have h := hφ.le_on_segment hr hys hseg
  have hry : φ r ≤ φ y := hmin y hys
  calc φ x ≤ max (φ r) (φ y) := h
    _ = φ y := max_eq_right hry

/-- A concave function on [lo, hi] with a maximum at r is nondecreasing on [lo, r]. -/
theorem concave_monotone_before_max {φ : ℝ → ℝ} {lo hi r : ℝ}
    (hφ : ConcaveOn ℝ (Icc lo hi) φ) (hr : r ∈ Icc lo hi) (hmax : ∀ z ∈ Icc lo hi, φ z ≤ φ r) :
    MonotoneOn φ (Icc lo r) := by
  intro x hx y hy hxy
  have hxs : x ∈ Icc lo hi := ⟨hx.1, le_trans hx.2 hr.2⟩
  have hys : y ∈ Icc lo hi := ⟨hy.1, le_trans hy.2 hr.2⟩
  have hseg : y ∈ segment ℝ x r := by
    rw [segment_eq_Icc (le_trans hxy hy.2)]
    exact ⟨hxy, hy.2⟩
  have h := hφ.ge_on_segment hxs hr hseg
  have hxr : φ x ≤ φ r := hmax x hxs
  calc φ x = min (φ x) (φ r) := (min_eq_left hxr).symm
    _ ≤ φ y := h

/-- Left flattening φ̃ = φ ∘ (max · r) is convex on [lo, hi]. -/
theorem flatten_left_convex {φ : ℝ → ℝ} {lo hi r : ℝ}
    (hφ : ConvexOn ℝ (Icc lo hi) φ) (hr : r ∈ Icc lo hi) (hmin : ∀ z ∈ Icc lo hi, φ r ≤ φ z) :
    ConvexOn ℝ (Icc lo hi) (φ ∘ fun x => max x r) := by
  have hmaxc : ConvexOn ℝ (Icc lo hi) (fun x : ℝ => max x r) :=
    (convexOn_id (convex_Icc lo hi)).sup (convexOn_const r (convex_Icc lo hi))
  have hsub : Icc r hi ⊆ Icc lo hi := fun z hz => ⟨le_trans hr.1 hz.1, hz.2⟩
  have hφ' : ConvexOn ℝ (Icc r hi) φ := hφ.subset hsub (convex_Icc r hi)
  have hmap : MapsTo (fun x : ℝ => max x r) (Icc lo hi) (Icc r hi) := by
    intro x hx
    exact ⟨le_max_right x r, max_le hx.2 hr.2⟩
  exact comp_convexOn_of_monotone hmaxc hφ' (convex_monotone_after_min hφ hr hmin) hmap

/-- Left flattening is nondecreasing on [lo, hi]. -/
theorem flatten_left_monotone {φ : ℝ → ℝ} {lo hi r : ℝ}
    (hφ : ConvexOn ℝ (Icc lo hi) φ) (hr : r ∈ Icc lo hi) (hmin : ∀ z ∈ Icc lo hi, φ r ≤ φ z) :
    MonotoneOn (φ ∘ fun x => max x r) (Icc lo hi) := by
  intro x hx y hy hxy
  have hm := convex_monotone_after_min hφ hr hmin
  exact hm ⟨le_max_right x r, max_le hx.2 hr.2⟩ ⟨le_max_right y r, max_le hy.2 hr.2⟩
    (max_le_max hxy le_rfl)

/-- Left flattening keeps the error band against a nondecreasing Γ. -/
theorem flatten_left_error {φ Γ : ℝ → ℝ} {lo hi r e : ℝ}
    (hr : r ∈ Icc lo hi) (hmin : ∀ z ∈ Icc lo hi, φ r ≤ φ z) (hΓ : MonotoneOn Γ (Icc lo hi))
    (herr : ∀ z ∈ Icc lo hi, |Γ z - φ z| ≤ e) :
    ∀ z ∈ Icc lo hi, |Γ z - (φ ∘ fun x => max x r) z| ≤ e := by
  intro z hz
  simp only [Function.comp]
  rcases le_total r z with h | h
  · rw [max_eq_left h]
    exact herr z hz
  · rw [max_eq_right h]
    have hrz := hmin z hz
    have hΓzr : Γ z ≤ Γ r := hΓ hz hr h
    have e1 := abs_le.mp (herr r hr)
    have e2 := abs_le.mp (herr z hz)
    rw [abs_le]
    constructor <;> linarith [e1.1, e1.2, e2.1, e2.2]

/-- Right flattening φ̃ = φ ∘ (min · r) is concave on [lo, hi]. -/
theorem flatten_right_concave {φ : ℝ → ℝ} {lo hi r : ℝ}
    (hφ : ConcaveOn ℝ (Icc lo hi) φ) (hr : r ∈ Icc lo hi) (hmax : ∀ z ∈ Icc lo hi, φ z ≤ φ r) :
    ConcaveOn ℝ (Icc lo hi) (φ ∘ fun x => min x r) := by
  have hminc : ConcaveOn ℝ (Icc lo hi) (fun x : ℝ => min x r) :=
    (concaveOn_id (convex_Icc lo hi)).inf (concaveOn_const r (convex_Icc lo hi))
  have hsub : Icc lo r ⊆ Icc lo hi := fun z hz => ⟨hz.1, le_trans hz.2 hr.2⟩
  have hφ' : ConcaveOn ℝ (Icc lo r) φ := hφ.subset hsub (convex_Icc lo r)
  have hmap : MapsTo (fun x : ℝ => min x r) (Icc lo hi) (Icc lo r) := by
    intro x hx
    exact ⟨le_min hx.1 hr.1, min_le_right x r⟩
  exact comp_concaveOn_of_monotone hminc hφ' (concave_monotone_before_max hφ hr hmax) hmap

/-- Right flattening keeps the error band against a nondecreasing Γ. -/
theorem flatten_right_error {φ Γ : ℝ → ℝ} {lo hi r e : ℝ}
    (hr : r ∈ Icc lo hi) (hmax : ∀ z ∈ Icc lo hi, φ z ≤ φ r) (hΓ : MonotoneOn Γ (Icc lo hi))
    (herr : ∀ z ∈ Icc lo hi, |Γ z - φ z| ≤ e) :
    ∀ z ∈ Icc lo hi, |Γ z - (φ ∘ fun x => min x r) z| ≤ e := by
  intro z hz
  simp only [Function.comp]
  rcases le_total z r with h | h
  · rw [min_eq_left h]
    exact herr z hz
  · rw [min_eq_right h]
    have hzr := hmax z hz
    have hΓrz : Γ r ≤ Γ z := hΓ hr hz h
    have e1 := abs_le.mp (herr r hr)
    have e2 := abs_le.mp (herr z hz)
    rw [abs_le]
    constructor <;> linarith [e1.1, e1.2, e2.1, e2.2]

/-- Error transfer through the probability-integral transform: F_n = Γ_n ∘ F. -/
theorem transfer_error {Γ ψ F : ℝ → ℝ} {e : ℝ}
    (hF : ∀ x, F x ∈ Icc (0 : ℝ) 1) (herr : ∀ s ∈ Icc (0 : ℝ) 1, |Γ s - ψ s| ≤ e) :
    ∀ x, |(Γ ∘ F) x - (ψ ∘ F) x| ≤ e := fun x => herr (F x) (hF x)

/-- DT7 Step E (validity direction of the LP): on a knot interval the data is a constant c and the
piecewise-linear fit is a convex combination of its end values, so the error inside is at most the
larger end error. Hence an LP-feasible knot vector gives a function whose sup error equals its knot
error, and the LP value is an upper bound for d(Γ_n, C_[a,b]). -/
theorem const_vs_linear_error (c ga gb t : ℝ) (ht0 : 0 ≤ t) (ht1 : t ≤ 1) :
    |c - ((1 - t) * ga + t * gb)| ≤ max |c - ga| |c - gb| := by
  have hsplit : c - ((1 - t) * ga + t * gb) = (1 - t) * (c - ga) + t * (c - gb) := by ring
  rw [hsplit]
  have h1 : |(1 - t) * (c - ga)| = (1 - t) * |c - ga| := by
    rw [abs_mul, abs_of_nonneg (by linarith)]
  have h2 : |t * (c - gb)| = t * |c - gb| := by
    rw [abs_mul, abs_of_nonneg ht0]
  calc |(1 - t) * (c - ga) + t * (c - gb)| ≤ |(1 - t) * (c - ga)| + |t * (c - gb)| := abs_add_le _ _
    _ = (1 - t) * |c - ga| + t * |c - gb| := by rw [h1, h2]
    _ ≤ (1 - t) * max |c - ga| |c - gb| + t * max |c - ga| |c - gb| :=
        add_le_add (mul_le_mul_of_nonneg_left (le_max_left _ _) (by linarith))
          (mul_le_mul_of_nonneg_left (le_max_right _ _) ht0)
    _ = max |c - ga| |c - gb| := by ring

/-- Size chain (GOAL Tier 3, gap G3): if the worst-case asymptotic size is 1 - F(2 q) (DT3 Thm 1'/Thm 4,
with F = F_M nondecreasing), q <= q_up (certified, DT8 Part A) and F(2 q_up) <= f_up (certified Arb
enclosure), then the size is at least 1 - f_up. -/
theorem size_chain (F : ℝ → ℝ) (hF : Monotone F) (q qup fup : ℝ) (hq : q ≤ qup)
    (hf : F (2 * qup) ≤ fup) : 1 - fup ≤ 1 - F (2 * q) := by
  have : F (2 * q) ≤ F (2 * qup) := hF (by linarith)
  linarith

/-- The DT8 numbers: q_.05 <= 0.55575 and F_M(1.112) <= 0.7160 (Arb), so the uniform-calibrated test has
worst-case asymptotic size >= 0.284 > 0.05 at the 5% level. -/
theorem size_chain_dt8 (F : ℝ → ℝ) (hF : Monotone F) (q : ℝ) (hq : q ≤ 0.55575)
    (hf : F 1.112 ≤ 0.7160) : (0.05 : ℝ) < 1 - F (2 * q) := by
  have h1 : F (2 * q) ≤ F 1.112 := hF (by linarith)
  linarith

/-! NONVAC: the hypotheses are satisfiable and the flattening really changes a non-monotone branch. -/
theorem nonvac_convex : ConvexOn ℝ (Icc (0 : ℝ) 1) (fun x : ℝ => (x - 1 / 2) ^ 2) := by
  refine ⟨convex_Icc 0 1, ?_⟩
  intro x _ y _ a b ha hb hab
  simp only [smul_eq_mul]
  have hb' : b = 1 - a := by linarith
  subst hb'
  nlinarith [mul_nonneg (mul_nonneg ha hb) (sq_nonneg (x - y))]

theorem nonvac_min : ∀ z ∈ Icc (0 : ℝ) 1, (fun x : ℝ => (x - 1 / 2) ^ 2) (1 / 2) ≤
    (fun x : ℝ => (x - 1 / 2) ^ 2) z := by
  intro z _
  show ((1 : ℝ) / 2 - 1 / 2) ^ 2 ≤ (z - 1 / 2) ^ 2
  nlinarith [sq_nonneg (z - 1 / 2)]

/-- The repaired branch is convex and nondecreasing, and it really differs from the original at 0. -/
example : ConvexOn ℝ (Icc (0 : ℝ) 1) ((fun x : ℝ => (x - 1 / 2) ^ 2) ∘ fun x => max x (1 / 2)) :=
  flatten_left_convex nonvac_convex ⟨by norm_num, by norm_num⟩ nonvac_min

example : MonotoneOn ((fun x : ℝ => (x - 1 / 2) ^ 2) ∘ fun x => max x (1 / 2)) (Icc (0 : ℝ) 1) :=
  flatten_left_monotone nonvac_convex ⟨by norm_num, by norm_num⟩ nonvac_min

example : ((fun x : ℝ => (x - 1 / 2) ^ 2) ∘ fun x => max x (1 / 2)) 0 = 0 ∧
    (fun x : ℝ => (x - 1 / 2) ^ 2) 0 = 1 / 4 := by norm_num

/-! Axioms audit (Step 10, gap G3). -/
#print axioms comp_convexOn_of_monotone
#print axioms comp_concaveOn_of_monotone
#print axioms flatten_left_convex
#print axioms flatten_left_monotone
#print axioms flatten_left_error
#print axioms flatten_right_concave
#print axioms flatten_right_error
#print axioms transfer_error
#print axioms const_vs_linear_error
#print axioms size_chain
#print axioms size_chain_dt8
