# DT15-A-PREREG — Step 16, work package A: C-edge∞ at the standard cutoffs (all three levels), v1

Status: **v1 FROZEN 2026-11-09 19:26Z** (`date -u` 19:25:26 at hashing). The v0 draft was written after `date -u` 19:21Z; v1 = v0 plus the fixes from one review round.

## 1. Claim
For each level α ∈ {0.10, 0.05, 0.01} with t′_α = q₋(α) − 10⁻⁴ (q₋ = lower end of the certified q̄ bracket, DT13 / DT8:
t′ = 0.4935, 0.5378, 0.6161 exactly, as decimal rationals):

  **g(u, t′_α) < g(1, t′_α) = 1 − F_M(2t′_α) for every u ∈ [½, 1)** (pointwise strict), hence (reflection)
  sup_{u∈(0,1]} g(u, t′_α) = g(1, t′_α) (EQUALITY; a strict sup over u < 1 is impossible since d_u → d₁ a.s.).
F_M is continuous (M = G(b₁) has a density; DT10 Lemma D / BP11), so P(M ≥ x) = 1 − F_M(x).

g(u, t) = P(d(H_u, K_½) ≥ t) is the Brownian-limit exceedance with the mode at quantile u (DT10-EDGE §0).
Consequence (endpoint enclosure; g non-increasing in t; q₋ ≤ q̄ ≤ q₊):
  1 − F_M(2q₊) ≤ sup_u g(u, q̄) ≤ 1 − F_M(2t′),
i.e. the mode-at-the-edge shape is the worst one in the Brownian limit at the real dip cutoffs, and through Theorem E this
gives the upper side of lim sup_n sup_F P(size). Tier after success: **conditional** on Lemmas S, F(d), Λ0, B′, G, D, J
(written; Lean for the deterministic parts is separate work) and the published inputs P1-P9, with
computational failure probability ≤ δ_total (§6). Not finite n.

## 2. Coverage ledger (u = 1 − ε², ε ∈ (0, 1/√2])
| piece | ε range | method | file |
|---|---|---|---|
| E0 | (0, 0.02] | analytic majorant H/ε² ≤ Φ(0.02) < 0, Arb | `code/dt16_neprime_arb.py` |
| E1 | [0.02, R_α] | NE′ boxes [e, 1.004e], Arb, contiguous from 0.02 | same |
| M | [ε_lo(α), 0.70710679] | paired two-bridge MC certificate, cells, sign test, IUT | `code/dt16_pair_cert.py` |
| refl | u < ½ | reflection lemma (DT10-EDGE) | — |
ε_lo(α) = 0.125 / 0.130 / 0.150 for α = 0.10 / 0.05 / 0.01 (chosen after DEV, disclosed §8). The level is covered iff
Φ < 0, R_α ≥ ε_lo(α) (E1 boxes contiguous; a closed-cell overlap), and the M decision is CERTIFIED. 0.70710679 > 1/√2.
Shared cell ends: every M cell is closed [e_i, e_{i+1}]; E1 boxes are closed; the union is closed.

## 3. The near-edge piece (E0, E1)
Bound (DT15-NEPRIME §4-§7, DT11 §5): g(1 − ε², t) − g(1, t) ≤ H(ε, t) := 2S(λt/ε) + Eψ(δX) − [F(x) − F(2t)],
ψ(w) ≤ min{F(x) − F(x − w), w[d(x − w)H_geom(w) + w d(x − w)d(x − 2w)]} (second branch only for x − 2w ≥ 0.5).
Validated inputs:
- F_M′ enclosure on cells of width 0.0025 on [0.9, 3.0) (`code/dt16_fm_density.py`, 3 renewal steps, Arb; deposit
  `deposits/dt16/fm_density.json`), widened outward by 10⁻⁹ relative (float conversions upstream).
- d(z): suffix maxima of Arb enclosures of v³f_E(v) on cells of width 0.002 on [0.5, 8], tail ≤ 4·10⁻⁴⁹ beyond 8.
- **H_geom upper step bound: CONFIRM Monte Carlo (§5b)**, exact law (P6 stick-breaking + Lemma H(a)), K = 160 faces,
  truncation term (2/3)^K·2/√(2π)·10⁶ ≤ 6·10⁻²³ added (an omitted face can only matter if its slope lies in the window
  (s_J, s_J + Δ_comp], Δ_comp ≤ 10⁶ else counted as A < r; omitted slopes are N(0, 1/L_i) given L, density ≤ √(L_i/2π);
  E Σ_{i>K} √L_i = 2(2/3)^K; on the complement Δ = Δ_comp and 1 − v_J ≥ computed, so conservative); 48 values r ∈ geomspace(10⁻⁵, 0.5); Clopper-Pearson one-sided at
  δ_H/48 each (Bonferroni), δ_H = 10⁻⁴; H_geom(w) ≤ bound at the smallest grid r ≥ w; 1 above 0.5.
- Boxes: H over [e_a, e_b] ≤ 2S(λt/e_b) + Σ_j T_b(r_j)[g(r_{j+1}) − g(r_j)] + T_b(r_J) − [CL(x_a) − CL(2t)],
  T_b(r) = exp(−2r(r + κ_b)), κ_b = (2 − λ)t/e_b, g = running max of the ψ upper bound at δ_b r (Stieltjes bound for a
  non-decreasing integrand under tail domination); λ ∈ {0.1, …, 1.0}, best per box.
- E0: Φ = 16.01 h⁴e^{−2h²}/(λt)² + L·H_geom(η)/(2a s₁) + L²ε₁²/(2a²s₁²) + e^{−B/ε₁²}/ε₁² − t·m, ε₁ = 0.02,
  a = (2 − λ)t, h = λt/ε₁ ≥ 2, L = d(2t − 2η), B = 2aηs₁ ≥ ε₁², m = min lower density on [2t, 2t/s₁];
  λ ∈ {0.5, 1, 1.5}, η ∈ {10⁻⁴, 3·10⁻⁴, 10⁻³, 3·10⁻³, 10⁻²}.
Deterministic given the H_geom file: one run, `python code/dt16_neprime_arb.py deposits/dt16/hgeom_confirm_r2026110848.json confirm`.

## 4. The middle piece (M)
Per path (two independent Brownian bridges b₁, b₂, DT10 Lemma S) and cell C = [e_i, e_{i+1}]:
Y_C = hi_C − lo₁ ∈ {−1, 0, 1}, E[Y_C] ≥ sup_{u∈C} P(d_u ≥ t′) − P(d_1 ≥ t′) (paired design; cell envelope):
- lo₁ = 1{G(b₁,PL)/2 − E₁′ ≥ t′ + 10⁻⁹} (G 2-Lipschitz in sup norm ⇒ G(b₁)/2 ≥ t′).
- hi_C = 0 iff the PL split path is certainly feasible (Lemma F(d), exact int64 predicates, unit 2⁻³², Λ evaluated at
  y′ < y; margins 4·10⁻⁹) at the **mixed corner** (a_lo = √(1 − e_i²) on side 1, e_{i+1} on side 2) at
  e = t′ − max(a_lo E₁′, e_{i+1}E₂′). Monotone in ε: aG₁ and aΛ₁(e/a) decrease, εG₂ and εΛ₂(e/ε) increase (Λ
  non-increasing, DT10 Λ0(ii)), so this implies d_u < t′ for every u in C (1-Lipschitz residual, Corollary S′).
- E_k′ = E_k(1 + 10⁻⁹) + 2⁻³²: E_k = sup|b_k − b_k,PL| drawn once, after stopping, from a separate RNG stream, from its
  exact CONDITIONAL law given the stopped skeleton, P(E ≤ x | S_stop) = K(x√m)^m (m = number of intervals; uniform
  dyadic refinement; observation-only stopping rule); 10⁻⁹ = slack for the float inverse; 2⁻³² > skeleton
  rounding |b − B/2³²| ≤ 2⁻³³, so E′ > the residual to the rounded skeleton strictly (strict d_u < t′).
- Λ ≥ 0 by definition (DT10 Λ0(i), max(0, ·)), needed for the corner monotonicity; the first condition
  2e ≥ max(a_lo G₁, e_{i+1}G₂) is checked before any Λ (it gives the Λ domain for every u in C); e ≤ 0 ⇒ not feasible.
- Refinement: b₁ only; dyadic midpoints N(mean of ends, h/4) (exact conditional law); decided from skeleton values only
  by the frozen rule `ambiguous` (window W(m) = (1 − 10⁻⁴)-quantile of E_m; cells with e_i < 0.2), m₀ = 2¹⁸,
  m_max = 2²²; a path still ambiguous at m_max is evaluated as is (conservative).
- Unresolved paths (|B| > 2³⁵) get hi = 1 in every cell, lo₁ = 0 (Y = +1); never dropped or redrawn.
- Cells: widths max(3·10⁻⁴, min(0.03, 0.08e²)) from ε_lo(α) to 0.70710679 (`cell_edges`).
Test per cell: exact sign test P_C = 2^{−K}Σ_{j≤n₊}C(K, j), K = n₊ + n₋ (Python integers), K = 0 ⇒ fail.
Decision per level: **CERTIFIED iff all shards complete and P_C < 10⁻³ for every cell of the level** (intersection-union
test, Berger 1982: false-certificate probability ≤ 10⁻³ per level, no Bonferroni). No partial release; a failed level is
"not certified", not "C-edge false".

## 5. Runs (fixed before any CONFIRM draw; one run each; crash ⇒ rerun the same shard command)
a. M: N = 33,000 paths, 11 shards × 3,000, **root 2026110847**:
   `python code/dt16_pair_cert.py 262144 4194304 33000 s 11 2026110847 confirm` (s = 0…10), then
   `python code/dt16_pair_cert_summary.py confirm 2026110847 1e-3`.
b. H_geom: N = 4·10⁷, **root 2026110848**: `python code/dt16_hgeom_cert.py 40000000 2026110848 confirm 1e-4`.
c. NE′: §3 command (after b).
RNG: numpy PCG64(SeedSequence([root, m₀, shard, 0/1])) for M; PCG64(SeedSequence([root, 16, 7])) for H_geom.

## 6. Error budget
δ_total per level ≤ 10⁻³ (M, IUT) + 10⁻⁴ (H_geom, all 48 r jointly) = 1.1·10⁻³; jointly over the three levels
≤ 3·10⁻³ + 10⁻⁴. Arb steps use no δ. Implementation-correctness is an assumption, checked by DEV tests (§8) and review.

## 7. Not covered
Finite n; Lean for E0/E1 assembly and for Lemmas G/D/B′ (separate step); the q̄ bracket width (enclosure only).

## 8. Development disclosure (all DEV roots banned for CONFIRM; CONVENTIONS L69)
- 2026110836 (low-t Δ sign, 66k pairs), 0837 (Lemma-B cert at low t: fails), 0839 (paired pilot m 2¹⁶/2²⁰),
  0840/0841 (slack Q, Lemma D checks), 0842 (H_geom continuum), 0843 (machine devtest 40 paths),
  0844 (machine v1 pilot 22k: sign-IUT reach 0.136/0.126/0.140), 0845 (H_geom 200k and 4·10⁶ dense: NE′ reach
  0.1359/0.1693/0.2185 with the dense grid), 0846 (machine v2 mixed corner, 4,400 paths: reach 0.119/0.119/0.140,
  P₊ = 0 for ε ≥ 0.11 at 10%).
- ε_lo(α) chosen after these: inside the DEV NE′ reach with margin (0.011 / 0.039 / 0.068) and where the v2 pilot has
  P₋ ≥ 0.002, P₊ ≈ 0. N chosen for power (≥ 40 expected nonzero Y at the lowest cell).

## 9. Freeze (sha256 prefixes at `date -u` 19:25:26)
`code/dt16_pair_cert.py` bc40d451; `code/dt16_pair_cert_summary.py` 65197a3e; `code/dt16_hgeom_cert.py` fe8430e5;
`code/dt16_neprime_arb.py` fa6706f6; `code/dt16_fm_density.py` 88a0ffd7; `deposits/dt16/fm_density.json` 2cb6cc02.
Smoke after the last code edit (19:25Z, DEV roots
0843/0845, outputs deleted): both programs run, summary prints COMPLETE/DECISION lines.
