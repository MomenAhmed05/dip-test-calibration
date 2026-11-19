# DT13-PREREG — Step 14: the finite-n size map. CONFIRM tier + EXACT n = 4.

Status: **v2.1 FROZEN 2026-11-07 ~22:57Z (`date -u` 22:56:17 at hashing), before any CONFIRM draw.** v2 was written after `date -u` 22:50:22; v2.1 adds five text fixes from a review round, with no design change. The v1 draft (22:44Z) was revised after an earlier review round. After the freeze nothing below may change; a change makes a new version and keeps this one.
DEV roots used for the design (CONVENTIONS CV-6): 2026110740 (n = 4 reduction), 0741 (n = 4 exact vs MC), 0742 (U1 pilot), 0743 (per-n search), 0744 (n = 4 case-split pilot), 0745 (K-step exact search seeds), 0746 (quick dry runs), 0747 (seeded errors), 0748 (Lemma N4 property test).

## 1. Objects (A-13)
- Statistic: D_alg = diptest 0.12.0 `dipstat(x, allow_zero=False)` = D_U (CV-1). Independent implementation: `lp_dip.py` (LP from the definition, all mode types).
- Null classes: U_ac ⊂ 𝒰_c (continuous unimodal CDFs) ⊂ U (CV-2, CV-4 v1.1). Every law drawn here is in U_ac.
- Calibration: q_n(1−α) = inf{t : H_n(t) ≥ 1 − α}, H_n = law of D under iid U(0, 1), simulated (CV-3). Rule: reject iff D ≥ q_n(1−α) (CV-4). The survival-p-value rule is NOT used as an equivalent (CV-4 v1.1).
- **Estimand: S_{n,α} = sup_{F ∈ U_ac} P_F(D_n ≥ q_n(1−α))** (the same bracket holds for 𝒰_c). At a named law: "rejection probability". The uniform anchor b_{n,α} = P_U(D_n ≥ q_n(1−α)) is estimated and reported; it equals α exactly at n = 4 (Part X). For n ≥ 5, b = α is NOT assumed: b is estimated by the Part W anchor f0 and reported with its enclosure whether or not that contains α. Ties in the reference at q̂ are recorded as a diagnostic only; they do not establish an atom of H_n at q.
- **Deliverable: a simultaneous confidence ENCLOSURE** L_{n,α} ≤ S_{n,α} ≤ U_{n,α}, plus the named-law map. A per-n maximiser is a WITNESS (a lower bound), never "the LF law".
- Two-step law W(w, ε) (CV-9): density (1 − w)/L on (0, L), w/ε on (L, 1), L = 1 − ε (unimodal for all w, ε). Ramp (CV-9): (1 − w)/L′ on (0, L′), then linear up to c_top on (L′, 1), L′ = 1 − ε_r, mass w on the ramp.

## 2. Scope of n
- The size map is made at a declared set of primary n. **Primary grid G = {4, 5, 6, 7, 8, 10, 12, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 300, 500, 750, 1000}** (21 values spanning 4 to 1000).
- On G: the two-sided enclosure L ≤ S ≤ U (Parts W, L, S, N4, X).
- Off G (every other integer n in [4, 1000]): the fixed-law map and the lower bound L only; the size there is enclosed only by [L, 1]. No interpolation claim. This limit is stated wherever the size map is reported.

## 3. Parts, laws, roots, code
Declared CONFIRM roots (each used once, banned after): **2026110750** (Part W reference), **2026110751** (Part W laws), **2026110752** (Part L laws), **2026110753** (Part S), **2026110754** (Part N4), **2026110755** (Part Q), **2026110756** (LP-audit selection).

### Part W — fixed witnesses and the uniform anchor (CONFIRM), `code/dt13_confirm.py`
- Reference: R = 1,000,000 uniform dips per n ∈ {4, …, 1000}, root 0750, spawn (0, n), one stream per n in batch order. Per (n, α): q̂ = X_(⌈(1−α)R⌉); bracket X_(l) ≤ q ≤ X_(u), binomial tails ≤ γ/(2·2991) each, γ = 0.002 (2991 = 997 × 3 pairs). Ties at q̂ are counted.
- Laws, N = 100,000 per (law, n), root 0751, spawn (1, f, n):
  - f0 U(0,1) (anchor b_{n,α}); f1 **W(0.1, 0.005) — PRIMARY fixed witness** (frozen at Step 3); f2 ramp (w = 0.1, ε_r = 0.01; frozen at Step 3); f3 W(0.4, 0.005) (secondary);
  - f4 **path witness W(w_n, 10⁻⁴), w_n = min(0.4, 1.5/√n)** — a frozen rule fitted on DEV 0743 (secondary; a fixed rule, not a per-n search).
- Per cell: counts at q̂, X_(l), X_(u); enclosure [CP lower at X_(u), CP upper at X_(l)], simultaneous over K_W = 5 × 997 × 3 = 14,955 cells at δ_W = 0.002 (two-sided Bonferroni). Class: ABOVE (lower > α) / BELOW (upper < α) / UNRESOLVED.

### Part L — per-n selected two-step laws (secondary; fresh CONFIRM re-estimates), `code/dt13_confirm.py`
- 63 streams = 21 n ∈ G × 3 selection levels; laws frozen in §10 (DEV 0743 best W(w, ε) per (n, α) on a 16 × 5 grid, `dev_search.json`@681ce72d; at n = 4 the Part X thin-top maximiser parameters at ε = 10⁻⁶). N = 100,000 each, root 0752, spawn (2, n, ia).
- Each stream is reported at all three α: **K_L = 63 × 3 = 189 cells, simultaneous at δ_L = 0.002** (two-sided Bonferroni over all 189). All 189 enter L.

### Part S — upper bound via Corollary U1′ (CONFIRM; "statistically certified"), `code/dt13_u1_cert.py`
- **Corollary U1′ (DT7-UNIF §1, PROVEN-written).** For every n, t and every partition {I_j} of [½, 1] into closed blocks: sup_{F∈𝒰_c} P_F(D_n ≥ t) ≤ max_j P(S̄_{n,I_j} ≥ t), where, for a uniform sample with empirical CDF Γ_n, S̄_{n,I} = max over the nonempty closed spacing pieces [a, b] = I ∩ [x_(k), x_(k+1)] (x_(0) = 0, x_(n+1) = 1) of d(Γ_n, C_[a,b]), C_[a,b] = {g continuous on [0, 1], convex on [0, b], concave on [a, 1]}, and d is the sup distance to both one-sided values of Γ_n. Proof chain: samplewise D(F_n) ≤ d(Γ_n, K_u) ≤ S̄_{n,I_j} for u = F(m) ∈ I_j (domination via φ̃∘F and H&H Thm 1; reflection for u < ½; blocks; pieces; exact knot LP).
- Notation: K_u = {g continuous on [0, 1]: convex on [0, u], concave on [u, 1]} (DT7 §0); C_[a,b] ⊆ K_u for every u ∈ [a, b].
- **Indicator contract (DT9 `dt9_block_cert.py`; frozen code).**
  - Knots t_0 < … < t_{K−1}: the set {0, x_(1), …, x_(n), a, b, 1} without duplicates; K = number of knots; i_a, i_b = indices of a and b. One-sided values lo_k = Γ_n(t_k−), hi_k = Γ_n(t_k).
  - LP (HiGHS, tolerance 10⁻¹⁰, fallbacks listed in the code): minimise e over (g_0, …, g_{K−1}, e) subject to hi_k − e ≤ g_k ≤ lo_k + e for all k; secant slopes nondecreasing on [t_0, t_{i_b}] (convex); nonincreasing on [t_{i_a}, t_{K−1}] (concave).
  - Certified upper bound from the RETURNED float vector g (no trust in the solver): ĝ = the piecewise-linear interpolant of g at the knots. band = max_k max(hi_k − g_k, g_k − lo_k) (Γ_n is constant between knots, so the interpolant's error is attained at knots). Slopes σ_k = (g_{k+1} − g_k)/(t_{k+1} − t_k); slope change J_k = σ_k − σ_{k−1} at knot k. Repair h = Σ_{convex-side k} δ_k·(t_k − s)⁺ − Σ_{concave-side k} δ′_k·(s − t_k)⁺ with δ_k = max(−J_k, 0) at the middle knots of the convex LP rows (1 ≤ k ≤ i_b − 1; these knots are ≤ a, since no sample point lies strictly inside a piece) and δ′_k = max(J_k, 0) at the middle knots of the concave LP rows (i_a + 1 ≤ k ≤ K − 2; these knots are ≥ b). For a = b the kink at a itself is free (C_[a,a] has no curvature condition there). Each term is 0 on [a, b], convex (resp. concave), and at most δ_k (resp. δ′_k) in absolute value on [0, 1]; ĝ + h is continuous, convex on [0, b] and concave on [a, 1], so ĝ + h ∈ C_[a,b], and sup|h| ≤ rep := Σδ_k + Σδ′_k. Hence d(Γ_n, C_[a,b]) ≤ band + rep.
  - Rounding allowance for evaluating band, σ, J and rep in IEEE-754 binary64 with round-to-nearest: rnd = 16u(Σ_k|σ_k| + K·rep + 1), u = 2⁻⁵³, plus SAFETY = 10⁻⁹ (covers threshold and band rounding, ~10⁻¹⁷). ub = band + rep + rnd + 10⁻⁹; a non-finite ub is an LP failure. Nodes of the branch and bound: ub of C_[a_i, b_j] bounds every piece inside (C_[a_i,b_j] ⊆ C_[a,b] for [a, b] ⊆ [a_i, b_j]).
  - Decision: the sample counts as NON-exceeding for a block only if every piece is shown to have ub ≤ X_(l) − 10⁻¹² (so S̄ < X_(l), strict, Corollary U1′); an unresolved node or a failed LP counts as an exceedance. Errors can only add exceedances, so the CP upper limit stays valid.
- Design: n ∈ G (21 values), N_S = 22,000 uniform samples each (11 shards × 10 chunks × 200), root 0753 (SeedSequence([root, n, shard, chunk])). Blocks [.5,.9], [.9,.97], [.97,.995], [.995,1]. Cutoffs: X_(l) of Part W per α (raw dip scale).
- U^S_{n,α} = max over the 4 blocks of the one-sided CP upper limit at confidence 1 − δ_S/252 (252 = 21 × 3 × 4), δ_S = 0.002.

### Part N4 — Lemma N4 upper bound at n = 4 (CONFIRM; statistically certified), `code/dt13_n4_bound.py`
- N = 20,000,000 sorted 4-tuples of uniforms, root 0754; 100 blocks of u on [½, 1]; block event E_I = E_L ∪ E_R (§4); one-sided CP upper limits at confidence 1 − δ_N4/300, δ_N4 = 0.001; U^{N4}_α = max over blocks. Cutoff: rational r_lo ≤ r_α (a larger event: conservative).

### Part X — exact n = 4 (EXACT tier, given the H&H reduction; no random numbers), `code/dt13_n4_exact.py`, checker `code/exact_poly_check.py`
- (a) H&H 1985 n = 4 reduction: T = max{1, s₂/max(s₁, s₃)}, D = T/(4(T+1)) **[LIT; S numerical check DEV 0740 (error ≤ 5.6e-17 on 20,000 samples); not re-proved here]**. (b) P_U(T ≥ r) = 2/((r+1)(r+2)) **for r > 1** (PROVEN: ∫₀^∞ e^{−y}(1 − e^{−y/r})² dy with iid exponential spacings; at r = 1 the left side is 1). Hence q_4(1−α) = r_α/(4(r_α+1)), (r_α+1)(r_α+2) = 2/α, r_.10 = 3; b_{4,α} = α exactly.
- Named two-step laws: P_F(T ≥ r) = Σ_k 24·c₁^k·c₂^{4−k}·vol(R_k(r)), R_k(r) = {0 ≤ x₁ ≤ … ≤ x₄ ≤ 1, x_k ≤ L ≤ x_{k+1} (k points below L), s₂ ≥ r s₁, s₂ ≥ r s₃}, c₁ = (1−w)/L, c₂ = w/ε; exact volumes (`exact_poly.py`: exact vertex enumeration + fan triangulation over tight-vertex faces; degenerate polytopes return 0). r_α irrational at α = .05, .01: rationals r_lo < r_α < r_hi (width < 10⁻³⁰), value interval [P(T ≥ r_hi), P(T ≥ r_lo)]; **the actual computed interval widths are reported** (not inferred from the r width).
- Checks: (i) the uniform as W(½, ½): at BOTH rational endpoints r_lo, r_hi the polytope value equals 2/((r+1)(r+2)) exactly, and the resulting probability interval contains α (checks regions and weights end to end; DEV run: true at all three levels, interval widths 1.1e-33 / 1.2e-32 / 3.5e-32); (ii) F-8 checker: Lasserre facet recursion (`exact_poly_check.py`), a different exact algorithm on the same regions (60/60 identical at DEV; re-run at the freeze); (iii) V2 (§5): the exact interval [p⁻, p⁺] of each named law must lie inside its fresh MC enclosure (checks regions, weights and the reduction end to end). Computed widths of all named-law intervals at the DEV run: ≤ 3.6e-32.
- Atom family F_p = pU + (1−p)δ₁ (atom row, CV-7): P = [p⁴ + 4p³(1−p)]α + 6p²(1−p)²/(1 + r_α) (S derivation; DEV 0741 MC within 1.21 SE at three points). Atom rows are never part of the U_ac size.
- **Thin-top limit value (NOT a supremum claim):** lim_{ε→0} P_{W(1−p, ε)}(T ≥ r_α) = f(p) = α(p⁴ + 4p³q + q⁴) + 6p²q²/(1 + r_α), q = 1 − p; max_p f = 3.1352% / 8.4593% / 13.6485% (α = .01/.05/.10). This is a LIMITING-FAMILY value and a lower bound on S_{4,α} (approached, not attained). That no finite-ε two-step law exceeds it is NOT claimed; the exact K-step search (K = 2, 3, 4; DEV 0745) found no higher value (A-3 near-miss record).
- The ramp at n = 4 is not computed exactly (polynomial density); MC only.

### Part Q — certified-MC bracket of q̄_.10 (+ replication of q̄_.05, q̄_.01), `code/qbar_cert.py` + `code/dt13_qbar_analysis.py`
- The DT8 Part A design unchanged (m = 2,000,000; N = 8,000,000; R = 12,000; 11 shards), root 0755, tag `q10cert`. Order-statistic brackets with ε = δ_Q/6 per side per α, δ_Q = 0.002. Asymptotic worst-case size 1 − F_M(2q̄_α) via the certified F_M enclosure (`fm_cert_h0.0002.npz`@943658cc; `fm_arb_low.npz`@0ab52bab below 0.9; both deterministic Arb enclosures with an F-8 checker, DT8/DT10).
- Replication: the fresh .05/.01 brackets must overlap DT8 Part A's; disagreement = escalate (a diagnostic, not a proof of correctness).

### LP audit (independent implementation)
- Every Part W law (f0-f4) and every Part L law, at n ∈ {4, 5, 6, 7, 8, 10, 20, 50, 100, 200} (Part L at those n in G): 100 random samples (root 0756) + every sample with |D − X| < 10⁻⁹ for any bracket value X. Rule: |D_alg − D_LP| ≤ 10⁻⁹ AND the two implementations agree on D ≥ X for every bracket value X. Any failure = audit FAIL.

## 4. Lemma N4 (n = 4 upper bound) — v2
**Statement.** Let F ∈ 𝒰_c with a mode m, u = F(m) ∈ [0, 1], and let a < b < c < d be the order statistics of 4 iid U(0,1); g₁ = b − a, g₂ = c − b, g₃ = d − c. Put A(v) = min(max(v − b, 0), g₂) and φ(v) = A(v)/g₁ + (g₂ − A(v))/g₃. **For every r > 1**, P_F(T ≥ r) ≤ P(φ(u) ≥ r). (φ(v) = g₂/g₃ for v ≤ b, = g₂/g₁ for v ≥ c, linear in between.) The statement is FALSE at r = 1 (F uniform, u = 1: P(T ≥ 1) = 1 but P(g₂ ≥ g₁) = ½; DEV 0748: 0.4991).
**Regularity.** X_(i) = H(U_(i)), H = F⁻¹ on (0, 1). On (−∞, m] F is convex and nondecreasing with F → 0 at −∞; if F were constant at a level c ∈ (0, u) on an interval [x₁, x₂] and later increased, then some chord to the left of x₁ has positive slope (F rises from 0 to c), the chord over [x₁, x₂] has slope 0, and a later chord has positive slope again: chord slopes of a convex function are nondecreasing, a contradiction. On [m, ∞) F is concave and nondecreasing; after a plateau its slope stays 0, so it cannot increase again below level 1. Hence F has no interior plateau and H is finite, continuous and strictly increasing on (0, 1), concave on (0, u] and convex on [u, 1). Unbounded densities, modal plateaus and u ∈ {0, 1} are allowed (only the relevant one-sided case is used; H(0), H(1) are never evaluated). σ_i = s_i/g_i > 0 are chord slopes.
**Proof.** Since r > 1, T ≥ r iff s₂ ≥ r s₁ and s₂ ≥ r s₃.
- u ≤ b: b, c, d ∈ [u, 1) where H is convex, so σ₂ ≤ σ₃; s₂ ≥ r s₃ gives σ₂g₂ ≥ rσ₃g₃ ≥ rσ₂g₃, so φ(u) = g₂/g₃ ≥ r.
- c ≤ u: a, b, c ∈ (0, u] where H is concave, so σ₁ ≥ σ₂; s₂ ≥ r s₁ gives φ(u) = g₂/g₁ ≥ r.
- b < u < c: concavity on [a, u] gives (H(u) − H(b))/(u − b) ≤ σ₁; convexity on [u, d] gives (H(c) − H(u))/(c − u) ≤ σ₃. So s₂ ≤ σ₁A + σ₃B, A = u − b, B = c − u. With λ = σ₃/σ₁ > 0 the two conditions give A + λB ≥ r g₁ and A/λ + B ≥ r g₃. If A ≥ r g₁ or B ≥ r g₃ then A/g₁ + B/g₃ ≥ r. Otherwise λ must lie in [(r g₁ − A)/B, A/(r g₃ − B)], nonempty iff (r g₁ − A)(r g₃ − B) ≤ AB iff φ(u) = A/g₁ + B/g₃ ≥ r.
- b = u or c = u: probability 0. ∎
**Reflection.** x ↦ −x maps F to F̃ ∈ 𝒰_c with ũ = 1 − u, swaps g₁ ↔ g₃ and keeps T; so S_{4,α} ≤ sup_{u∈[½,1]} P(φ(u) ≥ r_α).
**Blocks.** φ is constant–linear–constant in v, hence monotone on [L, R]; so ∪_{v∈[L,R]} {φ(v) ≥ r} = E_L ∪ E_R with E_v = {φ(v) ≥ r}. Therefore S_{4,α} ≤ max_j P(E_{L_j} ∪ E_{R_j}).
**Property test (DEV 0748, `code/dt13_n4_property.py`@9ade51e8).** 10⁷ checks of "T(H(U)) ≥ r ⇒ φ(u) ≥ r" over random concave/convex H (slopes over 12 decades), u ∈ {0, 1, random, near ½, near a sample point}, r ∈ {1 + 10⁻⁹, 1.5, 3, r_.05, r_.01}: 1,221,033 rejections, 0 violations; block containment 0 violations.
**Status.** PROVEN (written) for r > 1; a sample-wise bound, not expected to be tight (DEV 0744/0748: ~10.9/24.8/35.1% vs the thin-top limit 3.1/8.5/13.6%).

## 5. Claims (frozen wording), error budget, stop rule
- **E1 (size enclosure, primary).** With simultaneous MC confidence ≥ 1 − (γ + δ_W + δ_L + δ_S + δ_N4) = 0.991: for every n ∈ G and α ∈ {.01, .05, .10}, L_{n,α} ≤ S_{n,α} ≤ U_{n,α}; for every other integer n ∈ [4, 1000], S_{n,α} ≥ L_{n,α}. L = max of (the Part W CP lower limits of f1-f4, the Part L CP lower limits of the 189 cells at that n, and at n = 4 the LOWER endpoint p⁻ of the EXACT interval of each named n = 4 law, never a midpoint or a rounded value); U = U^S on G, and at n = 4 min(U^S, U^{N4}). Every number carries its tier label: "MC, statistically certified (simultaneous CP)" or "EXACT (given the H&H reduction)".
- **E2 (fixed-witness map).** Part W enclosures and classes for f1-f4 at every n and α (descriptive; no hypothesis "the witness over-rejects at every n").
- **E3 (n = 4, EXACT given the H&H reduction).** Named-law values with their computed interval widths; uniform b_{4,α} = α; the atom-family formula (atom row); the thin-top limit value (a lower bound, not a supremum). n = 5-8: no exact value (§8).
- **E4 (LF description, descriptive).** The Part L selections and the path rule with their CONFIRM values; never "the" LF law.
- **E5 (Part Q, separate family at δ_Q = 0.002).** q̄_.10 bracket and 1 − F_M(2q̄_.10); replication of q̄_.05, q̄_.01. A joint statement of E1 and E5 has confidence ≥ 0.989.
- **Validators (P-3).** V1: the uniform-anchor enclosure contains α in every (n, α) cell; any failure falls under the stop rule (no atom exception; reference ties are a diagnostic only). V2: at n = 4, for every law with an EXACT interval (f1, f3, f4, the three Part L laws), [p⁻, p⁺] ⊆ [MC lower, MC upper] (containment, not overlap). V6: L ≤ U on G. V8: k_{X_u} ≤ k_q̂ ≤ k_{X_l} in every cell. LP audit passes. Part Q replication overlaps DT8. Seeded errors (`code/dt13_seeded.py`, DEV 0747): 14/14 caught, clean pipeline passes all ten validators (`deposits/dt13/seeded_errors.json`@59ee41a4); this measures sensitivity to those errors only.
- **Stop rule.** Any validator failure stops the reading of E1-E5: no claim is made. There is no silent rerun with new seeds; a rerun is a new prereg version with new roots, and both results are reported. The stated confidence is unconditional; it is not reinterpreted as conditional on passing the validators.

## 6. S predictions (recorded before the run)
1. V1, V2, V6, V8, the audit and the Part Q replication pass.
2. W(0.1, 0.005) at 5%: 5.55% at n = 4 (exact), rising to ~25% at n = 1000; first ABOVE at 5% and 10% by n ≈ 10-20; at 1% by n ≈ 20-40.
3. Path witness: L ≈ 3.0/8.4/13.6% at n = 4 and ≈ 13/29/40% at n = 1000 (1%/5%/10%).
4. Part L fresh values 0-2 points below their DEV values (winner's curse).
5. U^S ≈ 14/29/39% at n = 1000, ≈ 19/39/52% at n = 50, above 25/44/56% at n ≤ 20; U^{N4} ≈ 11/25/35%.
6. "Informative" := U − L ≤ 10 points. Expected informative only at n ≥ 500 (5%, 10%) and n ≥ 200 (1%).
7. q̄_.10 ≈ 0.50; asymptotic worst-case size at 10% ≈ 41-47%; replication brackets overlap DT8.

## 7. What each outcome means
- Witness lower > α: nominal violation at that cell (not global LF). Witness above the anchor: the uniform is not least favourable there.
- An enclosure straddling α: inconclusive about validity (never evidence of validity). U ≤ α would show validity (not expected).
- A narrow enclosure with unstable winning parameters: size determined, LF configuration not identified.
- A validator failure: §5 stop rule.

## 8. Not done, and why (feasibility record; A-9 check)
- **Exact n = 5-8 is NOT pursued under this budget (scope reduction, disclosed).** Reason: we have not established a usable polyhedral decomposition of {D ≥ c} for n = 5-8; the decoupled pinned formulation we tested is not equivalent to D (it differs by up to 0.035, DEV 0740), and the one integration strategy costed (certified boxes) needs ~10⁶ boxes for ±0.1 point at n = 5 (uniform only, n − 2 dimensions) and ≥ 10⁹ for n ≥ 6 or for two-step laws (n dimensions, no affine reduction). This is "not pursued under this budget", not "shown infeasible": a cheaper symbolic decomposition may exist. Rows n = 5-8: MC with exact binomial bounds + LP audit (every random and near-threshold sample, decision agreement).
- **A tight finite-n upper bound at small n.** U1′ and Lemma N4 are sample-wise bounds; at small n they are 2.5-10× the best law found (DEV 0742/0744). A tight bound needs a non-sample-wise argument: OPEN.

## 10. Frozen Part L table (from `dev_search.json`@681ce72d; n = 4 from Part X)
n: (w, ε) for α = .01 | .05 | .10
- 4: (0.473789162301, 1e-6) | (0.433348270040, 1e-6) | (0.396743906214, 1e-6)
- 5: (0.4, 1e-6) | (0.4, 0.02) | (0.4, 0.02);  6: (0.45, 0.001) | (0.4, 0.001) | (0.4, 0.001);  7: (0.4, 0.02) | (0.4, 0.02) | (0.35, 1e-6);  8: (0.35, 1e-6) ×3
- 10: (0.4, 0.005) | (0.3, 1e-6) | (0.3, 0.001);  12: (0.35, 0.001) | (0.35, 0.005) | (0.3, 1e-4);  15: (0.35, 0.02) | (0.3, 0.001) | (0.3, 1e-4)
- 20: (0.3, 0.001) | (0.3, 1e-4) | (0.25, 1e-6);  25: (0.3, 1e-4) | (0.25, 1e-6) | (0.25, 1e-6);  30: (0.25, 0.001) | (0.25, 0.001) | (0.2, 1e-6)
- 40: (0.25, 0.001) | (0.2, 0.005) | (0.2, 0.005);  50: (0.2, 0.005) ×3;  75: (0.2, 1e-6) | (0.15, 1e-6) | (0.15, 1e-6);  100: (0.15, 0.001) ×3
- 150: (0.15, 1e-6) | (0.1, 0.001) | (0.1, 1e-6);  200: (0.1, 1e-4) ×3;  300: (0.1, 0.001) | (0.1, 0.001) | (0.07, 1e-4)
- 500: (0.07, 1e-4) | (0.07, 0.001) | (0.07, 0.001);  750: (0.07, 1e-4) | (0.05, 0.001) | (0.05, 0.001);  1000: (0.05, 1e-6) | (0.05, 1e-4) | (0.05, 1e-4)

## 11. Code, data and run order
Frozen sha256 prefixes (22:56Z): dt13_confirm.py cca8dd45; dt13_u1_cert.py 4c9fb8ba; dt13_n4_bound.py e26de0f0; dt13_n4_property.py 9ade51e8; dt13_analysis.py ca3592cc; dt13_qbar_analysis.py 4670b37a; qbar_cert.py 9506890b; dt9_block_cert.py cd8eaa67; sup_pinned.py 2f2441b1; pinned_dip.py 791dff28; dt0c_screen.py 7fe20a4f; lp_dip.py 6a00c686; exact_poly.py 4878f18a; exact_poly_check.py 8fb43481; dt13_n4_exact.py de01a942; dt13_seeded.py b94dd413; dt8_cert_analysis.py bf02a7ff; dev_search.json 681ce72d; n4_exact.json a74d8b9f; seeded_errors.json 59ee41a4; fm_cert_h0.0002.npz 943658cc; fm_arb_low.npz 0ab52bab; CONVENTIONS.md ebf3488c; DT7-UNIF.md 13587899. Freeze re-runs: Lasserre checker 60/60 PASS; seeded errors 14/14, clean passes.
Code: dt13_confirm.py, dt13_u1_cert.py, dt13_n4_bound.py, dt13_n4_property.py, dt13_analysis.py, dt13_qbar_analysis.py, qbar_cert.py, dt9_block_cert.py, sup_pinned.py, pinned_dip.py, dt0c_screen.py, lp_dip.py, exact_poly.py, exact_poly_check.py, dt13_n4_exact.py, dt13_seeded.py, dt8_cert_analysis.py. Data: dev_search.json, n4_exact.json, seeded_errors.json, fm_cert_h0.0002.npz, fm_arb_low.npz.
Run order: (1) Part X `python code/dt13_n4_exact.py` + `python code/exact_poly_check.py` (deterministic); (2) Parts W/L + audit `python code/dt13_confirm.py`; (3) Part S `python code/dt13_u1_cert.py n 10 s 11 2026110753 cert 200` for every n ∈ G and s = 0..10; Part N4 `python code/dt13_n4_bound.py 2026110754 20000000`; (4) Part Q `python code/qbar_cert.py 2000000 8000000 12000 s 11 2026110755 q10cert` for s = 0..10; (5) `python code/dt13_analysis.py 2026110753 cert 2026110754`; `python code/dt13_qbar_analysis.py 2026110755 q10cert`.

## 12. In simple words
- The question: the standard dip test says "5%". How large can its real error rate be, for sample sizes from 4 to 1000, at 1%, 5% and 10%?
- We cannot test every one-peak law. So at 21 chosen sample sizes we give two numbers. The low number comes from real laws, tested on new random data. The high number comes from a theorem that holds for every law at once. The true worst case lies between them.
- At every other sample size from 4 to 1000 we give the low number only.
- At n = 4 we also compute exact values with fractions, with no random numbers.
- We write all choices down now, before we draw the data, so we cannot pick them after we see the results.
