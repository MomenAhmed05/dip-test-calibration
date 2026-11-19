# DT10-PREREG — Step 11, task 2: certificate of C-edge∞ on the middle range of mode positions

Status: **FROZEN** 2026-11-07 ~12:07Z (`date -u`; draft ~12:00Z, numerics revision ~12:05Z), after two review rounds.

## Claim to certify
For α ∈ {0.05, 0.01}, c\* = ½F_M⁻¹(1 − α), and every ε-interval I_i of the fixed grid below that passes:

  sup_{ε ∈ I_i} g(1 − ε², c\*) < α = g(1, c\*),

with confidence ≥ 0.99 jointly over all 72 interval × level claims. Here g(u, t) is the Brownian-limit exceedance probability with the mode at quantile u (DT10-EDGE §0). Reflection covers u < ½.

## Why this is valid (all written in DT10-EDGE §6)
1. Lemma B: g(1 − ε², t) ≤ B̃(ε, t).
2. B̃ is non-increasing in t, and c\* ≥ t_lo (rigorous Arb bracket, DT8-RESULT: t_lo = 0.7305 for α = 0.05, 0.8585 for α = 0.01). So B̃(ε, c\*) ≤ B̃(ε, t_lo).
3. g(1, c\*) = 1 − F_M(2c\*) = α exactly, by the definition of c\* (F_M is continuous, Lemma D).
4. For ε ∈ I_i = [e_i, e_{i+1}], B̃(ε, t_lo) ≤ E[V_i] (pathwise), with
   - V_i = 1 if e_{i+1}(G_m + 2E) ≥ 2t_lo;
   - otherwise V_i = 1 − F_lo((2t_lo − e_{i+1}Λ⁺)/a_i), with a_i = √(1 − e_i²), Λ⁺ = Λ_{b_m}(t_lo/e_{i+1} − E) + E, argument reduced by 1e-9;
   - b_m = piecewise-linear Brownian bridge on m = 131072 equal subintervals; E = ‖b − b_m‖ drawn exactly and independently (m^{−1/2}·max of m Kolmogorov suprema); G_m = G(b_m); |G(b) − G(b_m)| ≤ 2E; Λ0(ii) monotone in y and Λ0(iv) perturbation;
   - F_lo = cumulative max of the Arb lower enclosures of F_M at the grid points of `deposits/dt10/fm_arb_low.npz` (floor lookup; 0 below 0.5). The table on [1.0, 8.5] is identical to the DT8 table `fm_arb_wide.npz`; [0.5, 1.0) is new (step 0.002; FE_ball remainder valid for z ≥ 0.3).
   - V_i is clipped to [0, 1]; a non-finite value is set to 1.
   - **Numerics (fixed before the freeze).** The skeleton is rounded to integers B_k = rint(2^40 b_k); the bridge used is B/2^40, and E is replaced by E′ = E + 2^−40. All hull, tangent and branch decisions for G and Λ use exact int64 arithmetic (`code/dt10_exact.py`; bounds |B| ≤ 2^43, 2y ≤ 16, products < 2^63). The level y is replaced by y′ = (floor(2^41 y₀) − 2)/2^41 < y₀, which is safe because Λ is non-increasing in y. Values are evaluated from exact integers (float error ≤ ~1e-14) and padded by +1e-12 (G and Λ). The trigger fires at ≥ 2t − 1e-12. F_lo values are lowered by 1e-12 (float(arb.lower) rounding). The CDF argument is lowered by 1e-9. A path with |B| > 2^43 or Y′ out of range gets V = 1 in every cell (counted as `unresolved_paths`); paths are never dropped or redrawn. The final comparison is upper bound + 1e-12 < α.
   - DEV check (root 2026110736): exact kernels = brute-force rational computation in 300 random cases (m ≤ 64; max |diff| 1.1e-16); at m = 131072 they agree with the old float kernels within 8e-13 (20 paths).
5. Empirical-Bernstein bound (Maurer & Pontil 2009, Thm 4) on E[V_i], d = 0.01/72 per claim (Bonferroni; shared paths across claims are allowed; replications i.i.d.).

## Design (fixed before any CONFIRM draw)
- Grid: K = 36 intervals, geometric from e_0 = 0.12 to e_36 = 0.70710679 (> 1/√2). u-range covered if all pass: [0.5, 0.9856].
- Levels: (α, t_lo) = (0.05, 0.7305), (0.01, 0.8585).
- n = 110,000 bridges (11 shards × 10,000), fixed; ONE look; no stopping rule.
- RNG: PCG64(SeedSequence([root, m, shard])), **CONFIRM root 2026110734**, shards 0–10.
- Concentration: s² is the unbiased sample variance (divisor n − 1) of the per-path V_i over all n = 110,000 paths pooled from the 11 shards; L = ln(2·72/0.01).
- Decision: interval i certified at level α iff upper bound + 1e-12 < α. The contiguous block is taken downward from the interval that contains 1/√2; it is empty if that interval fails. No bridging of gaps, no extrapolation below 0.12. No root replacement, path rejection or reruns other than crash-resume with the same seeds. Report: every certified interval, and the largest contiguous certified block ending at 1/√2. Failed intervals are "not certified", not "C-edge false".
- Commands: `python code/dt10_cert.py 131072 110000 s 11 2026110734 cert` for s = 0…10, then `python code/dt10_cert_analysis.py 2026110734 cert 131072`. On a crash: rerun the same shard command (same seeds).

## DEV disclosed (not CONFIRM)
- Root 2026110731: Brownian pair DEV (16k, 2,000 paths). Root 2026110732: float-F_M B(ε) scan (m 4096, 20,000 paths). Root 2026110733: V_i with m = 131072, 4,000 paths, ε ∈ [0.2, 1/√2] (dev2 with the old table; dev3 with the 0.5 table): min mean margin 0.0025 (α = 0.05, ε ≈ 0.2), 0.0010 (α = 0.01); sd ≤ 0.13. Root 2026110735: pipeline smoke (220 paths).
- The lower end e_0 = 0.12 was chosen after DEV; intervals below ≈ 0.15 are expected to fail at α = 0.05. This is disclosed; every interval is a separate pre-registered claim.

## What this does NOT give
- The band u ∈ (1 − e_0², 1): only Theorem N (strict local maximum, non-explicit ε₀(t)). C-edge∞ stays OPEN there until a quantitative Lemma D / Lemma Λ closes it.
- Finite n: this is the Brownian limit only (Theorem E links it to lim sup_F).

## Freeze (sha256, first 8)
- code/dt10_cert.py 45d9970f; code/dt10_exact.py dae342bb; code/dt10_cert_analysis.py 1b9f7280; code/dt10_edge_mc.py 9efb4043; code/fm_cert_check2.py b621d4bd; code/fm_cert_check.py 17255eda; deposits/dt10/fm_arb_low.npz 0ab52bab.
