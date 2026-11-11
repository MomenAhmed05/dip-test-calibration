/-
DipCore.lean — machine checks of the inequality steps in work/DT3-THEOREM.md (v0.5).
Step 8 (DT6) first item, started 2026-11-06 while the DT3c CONFIRM run is in progress.

1. `pinned_reduction`   : Lemma 5 / Lemma 7 lower bound, near-minimiser form (repair).
2. `ext_error_left`     : Lemma 7 step 5, extension to the left of block j stays within r*E of the data.
3. `ext_error_right`    : Lemma 7 step 5, mirror (right of block j, up to the turn point).
4. `seam_select`        : Lemma 7 step 6, the modal piece wins the max at the turn point.
5. `tail_split`         : Theorem 1' step 4 (repair), as a measure inequality.
6. `kuiper_eq_two_dist` : range(h) = 2 * inf_c sup |h - c| (used for the Kuiper = 2 x dip identity), finite form.
-/
import Mathlib

open MeasureTheory

/-- **Pinned reduction (Lemma 5, near-minimiser form).**
`e g` is the sup-error of a candidate fit `g` with turn point `m g`; `D` is a lower bound reached
arbitrarily closely (`hD`: every level above `D` is beaten by some fit).
If every fit with turn point `≤ t` has error `≥ κ - ε`, every fit with turn point `> t` has error `≥ V`,
`D ≤ ε` and `ε < κ - ε`, then `V ≤ D`. -/
theorem pinned_reduction {G : Type*} (e : G → ℝ) (m : G → ℝ) (t κ ε V D : ℝ)
    (hD : ∀ y, D < y → ∃ g, e g < y)
    (hleft : ∀ g, m g ≤ t → κ - ε ≤ e g)
    (hright : ∀ g, t < m g → V ≤ e g)
    (hDε : D ≤ ε) (hgap : ε < κ - ε) : V ≤ D := by
  by_contra hlt
  rw [not_le] at hlt
  -- pick a fit below min V (κ - ε)
  obtain ⟨g, hg⟩ := hD (min V (κ - ε)) (lt_min hlt (lt_of_le_of_lt hDε hgap))
  rcases le_or_gt (m g) t with h | h
  · have := hleft g h
    have : e g < κ - ε := lt_of_lt_of_le hg (min_le_right _ _)
    linarith
  · have := hright g h
    have : e g < V := lt_of_lt_of_le hg (min_le_left _ _)
    linarith

/-- **Extension error, left side (Lemma 7 step 5).**
Block `j` starts at `a`; the point `x ≤ a` lies to its left. `F` grows at most at rate `c - Δ` on `[x, a]`,
the fit `g` is within `E` of the smoothed target `h` at `a`, and `h` is `K`-Lipschitz.
The left extension `A x = F a + r g a - (c - r K) (a - x)` exceeds the data `F x + r h x`
by at most `r E - (Δ - 2 r K) (a - x)`. -/
theorem ext_error_left (Fx Fa ga ha hx r c Δ K E x a : ℝ) (hr : 0 ≤ r)
    (hF : Fa - Fx ≤ (c - Δ) * (a - x))
    (hg : |ga - ha| ≤ E)
    (hh : |ha - hx| ≤ K * (a - x)) :
    (Fa + r * ga - (c - r * K) * (a - x)) - (Fx + r * hx)
      ≤ r * E - (Δ - 2 * r * K) * (a - x) := by
  have h1 : ga - ha ≤ E := (abs_le.mp hg).2
  have h2 : ha - hx ≤ K * (a - x) := (abs_le.mp hh).2
  have h3 : r * (ga - ha) ≤ r * E := mul_le_mul_of_nonneg_left h1 hr
  have h4 : r * (ha - hx) ≤ r * (K * (a - x)) := mul_le_mul_of_nonneg_left h2 hr
  nlinarith [h3, h4, hF]

/-- **Extension error, right side (Lemma 7 step 5, mirror).** Block `j` ends at `b`; `x ≥ b`.
`F` grows at least at rate `c + Δ` on `[b, x]`. The right extension
`A x = F b + r g b + (c + r K) (x - b)` exceeds the data by at most `r E - (Δ - 2 r K) (x - b)`. -/
theorem ext_error_right (Fx Fb gb hb hx r c Δ K E x b : ℝ) (hr : 0 ≤ r)
    (hF : (c + Δ) * (x - b) ≤ Fx - Fb)
    (hg : |gb - hb| ≤ E)
    (hh : |hb - hx| ≤ K * (x - b)) :
    (Fb + r * gb + (c + r * K) * (x - b)) - (Fx + r * hx)
      ≤ r * E - (Δ - 2 * r * K) * (x - b) := by
  have h1 : gb - hb ≤ E := (abs_le.mp hg).2
  have h2 : hb - hx ≤ K * (x - b) := (abs_le.mp hh).2
  have h3 : r * (gb - hb) ≤ r * E := mul_le_mul_of_nonneg_left h1 hr
  have h4 : r * (hb - hx) ≤ r * (K * (x - b)) := mul_le_mul_of_nonneg_left h2 hr
  nlinarith [h3, h4, hF]

/-- **Seam selection (Lemma 7 step 6).** At the turn point `μ`, a non-modal extension is at most
`data + r E - (Δ - 2 r K) η` (with `μ - a_m ≥ η`, by `ext_error_*`), and the modal piece is at least
`data - r E`. If `2 r E < (Δ - 2 r K) η`, the modal piece is strictly larger. -/
theorem seam_select (Aj Am data r E Δ K η : ℝ)
    (hAj : Aj - data ≤ r * E - (Δ - 2 * r * K) * η)
    (hAm : data - r * E ≤ Am)
    (hsmall : 2 * r * E < (Δ - 2 * r * K) * η) : Aj < Am := by
  linarith

/-- **Tail split (Theorem 1' step 4 repair).** For a measure `P`, events
`A = {D ≥ q}`, `B = {V > q}`, `N = {V = q}` with `A ⊆ B ∪ N` (because `D ≤ V`) and `P N = 0`:
`P A + P (B \ A) ≤ P B`. With `P A ≥ α` this gives `P B ≥ α + P (B \ A)`. -/
theorem tail_split {Ω : Type*} [MeasurableSpace Ω] (P : Measure Ω) (A B N : Set Ω)
    (hA : MeasurableSet A) (hsub : A ⊆ B ∪ N) (hN : P N = 0) :
    P A + P (B \ A) ≤ P B := by
  have hAB : P A ≤ P (B ∩ A) := by
    calc P A ≤ P ((B ∩ A) ∪ N) := by
            apply measure_mono
            intro ω hω
            rcases hsub hω with h | h
            · exact Or.inl ⟨h, hω⟩
            · exact Or.inr h
      _ ≤ P (B ∩ A) + P N := measure_union_le _ _
      _ = P (B ∩ A) := by rw [hN, add_zero]
  calc P A + P (B \ A) ≤ P (B ∩ A) + P (B \ A) := by gcongr
    _ = P B := measure_inter_add_sdiff B hA

/-- **Range = 2 × best constant fit**, for a function on a finite nonempty set
(the Kuiper = 2 × dip identity uses this pointwise in the class, see DT5-DESIGN §3b).
For every `c`, `sup h - inf h ≤ 2 * max |h - c|`, with equality at the midpoint `c`. -/
theorem range_le_two_dist {ι : Type*} [Fintype ι] [Nonempty ι] (h : ι → ℝ) (c : ℝ) :
    (Finset.univ.sup' Finset.univ_nonempty h) - (Finset.univ.inf' Finset.univ_nonempty h)
      ≤ 2 * Finset.univ.sup' Finset.univ_nonempty (fun i => |h i - c|) := by
  obtain ⟨i, -, hi⟩ := Finset.exists_mem_eq_sup' Finset.univ_nonempty h
  obtain ⟨j, -, hj⟩ := Finset.exists_mem_eq_inf' Finset.univ_nonempty h
  have hi' : |h i - c| ≤ Finset.univ.sup' Finset.univ_nonempty (fun i => |h i - c|) :=
    Finset.le_sup' (fun i => |h i - c|) (Finset.mem_univ i)
  have hj' : |h j - c| ≤ Finset.univ.sup' Finset.univ_nonempty (fun i => |h i - c|) :=
    Finset.le_sup' (fun i => |h i - c|) (Finset.mem_univ j)
  rw [hi, hj]
  have a1 := le_abs_self (h i - c)
  have a2 := neg_abs_le (h j - c)
  linarith

/-! ### NONVAC: the hypotheses are satisfiable (concrete instances) -/

/-- `pinned_reduction` hypotheses hold for a two-fit family: a wrong-mode fit (`false`, turn point -1)
with error 1 ≥ κ - ε = 0.8, and a right-mode fit (`true`, turn point 1) with error 0.2 = V = D = ε;
t = 0, κ = 1. (A first attempt with D = 0.1 was rejected by Lean: D must be the infimum, 0.2.) -/
example : ∃ (e m : Bool → ℝ) (t κ ε V D : ℝ),
    (∀ y, D < y → ∃ g, e g < y) ∧ (∀ g, m g ≤ t → κ - ε ≤ e g) ∧
    (∀ g, t < m g → V ≤ e g) ∧ D ≤ ε ∧ ε < κ - ε ∧ V = 1/5 := by
  refine ⟨fun b => if b then 1/5 else 1, fun b => if b then 1 else -1, 0, 1, 1/5, 1/5, 1/5, ?_, ?_, ?_,
    le_refl _, by norm_num, rfl⟩
  · intro y hy
    exact ⟨true, by simpa using hy⟩
  · intro g hg
    cases g
    · norm_num
    · norm_num at hg
  · intro g hg
    cases g
    · norm_num at hg
    · norm_num

/-- `seam_select` hypotheses are satisfiable (r small relative to the kink). -/
example : ∃ Aj Am data r E Δ K η : ℝ, Aj - data ≤ r * E - (Δ - 2 * r * K) * η ∧
    data - r * E ≤ Am ∧ 2 * r * E < (Δ - 2 * r * K) * η :=
  ⟨-1, 0, 0, 1/100, 1, 1, 1, 1/10, by norm_num, by norm_num, by norm_num⟩

theorem two_dist_mid_le_range {ι : Type*} [Fintype ι] [Nonempty ι] (h : ι → ℝ) :
    2 * Finset.univ.sup' Finset.univ_nonempty
        (fun i => |h i - (Finset.univ.sup' Finset.univ_nonempty h
                          + Finset.univ.inf' Finset.univ_nonempty h) / 2|)
      ≤ (Finset.univ.sup' Finset.univ_nonempty h) - (Finset.univ.inf' Finset.univ_nonempty h) := by
  set S := Finset.univ.sup' Finset.univ_nonempty h
  set I := Finset.univ.inf' Finset.univ_nonempty h
  have key : ∀ i ∈ (Finset.univ : Finset ι), |h i - (S + I) / 2| ≤ (S - I) / 2 := by
    intro i _
    have hs : h i ≤ S := Finset.le_sup' h (Finset.mem_univ i)
    have hl : I ≤ h i := Finset.inf'_le h (Finset.mem_univ i)
    rw [abs_le]
    constructor <;> linarith
  have := Finset.sup'_le Finset.univ_nonempty _ key
  linarith

/-! ### Axioms audit -/
#print axioms pinned_reduction
#print axioms ext_error_left
#print axioms ext_error_right
#print axioms seam_select
#print axioms tail_split
#print axioms range_le_two_dist
#print axioms two_dist_mid_le_range
