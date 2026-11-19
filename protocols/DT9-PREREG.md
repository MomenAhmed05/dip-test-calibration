# DT9-PREREG — Step 11, task 1: larger certificate of c = 0.7315 at n = 800 and 1600, and a certified replay of DT8

Status: **FROZEN** 2026-11-07 ~01:57Z (`date -u`). Drafted ~01:40Z, revised ~01:55Z after two review rounds. Purpose: (1) a larger run to certify 0.7315 at n = 800 and 1600; (2) a certified replay of DT8.

## Why
DT8 certified c = 0.7315 (α = 0.05) for n ≤ 400, but only c = 0.74 at n = 800 (N = 30k) and n = 1600 (N = 15k). "Not certified" is not "invalid". The DT8 edge-block rate at 0.7315 was 4.86% (n = 800) and 4.73% (n = 1600), too close to 5% for those N.
A review also found a numerical gap in DT8: its LP values had no repair term, and the 1e-9 margin was not proved to cover shape violations. Part R below replays DT8 with certified numerics.

## Disclosure: what was seen before freezing (DEV / planning only)
- DT8 CONFIRM data (root 2026110711) were used to plan N and the cutoffs. Part A claims rest on the fresh root only.
- DEV root 2026110720 (`code/dt9_dev_check.py`, `deposits/dt9/dev_check.json`):
  - decision vs brute force: 4,100 checks, 6 disagreements. All are the safe kind: the cutoff lies 2–5e-11 below a single piece's own bound, and a node bound ≤ cutoff proves "no".
  - exact rational repair test: 80 perturbed LP solutions, all with violations; 0 shape errors and 0 bound violations.
  - injected LP failures at rates 0.1 / 0.5 / 1.0: 0 exceedances lost at the fixed cutoff. Three safe-kind disagreements at the 1e-9-relative cutoff (gaps 3–4e-11).
  - a cutoff equal to the value is never counted; near-coincident points (1e-13 apart) and a point exactly on a block boundary run without error.
- DEV root 2026110721 (driver smoke test, 11 × 100 samples per n): 0.21 s/sample at n = 800 and 0.57 s/sample at n = 1600 (11 processes). Repair term up to 2.4e-6 raw (≈ 7e-5 on the √n scale).
- DEV root 2026110723: full pipeline test (orchestrator + driver + analysis) with 1-chunk looks; the stop rule was applied correctly at look 1.
- DT8 sensitivity audit: if every DT8 sample within 1e-4 (√n scale) below a cutoff were an exceedance, all 12 DT8 cells still pass. This is sensitivity evidence only; Part R replaces it.

## Common computation (one-sided numerics)
- Per sample and block: a **conservative indicator** Z ≥ 1{√n S̄ > c} by decision branch and bound (`code/dt9_block_cert.py`). Z = 0 only if a certified upper bound proves every piece ≤ c/√n; everything else (exceedance or unresolved) is Z = 1. Then {Σ Z ≤ b} ⊆ {Σ Y ≤ b} pathwise, so the binomial calibration of the true indicators applies.
- Node bound: an upper bound on d(Γ_n, C_[a_i, b_j]), which is ≥ every piece inside, because C_[A,B] ⊆ C_[a,b] for [a,b] ⊆ [A,B].
- **Certified primal upper bound** of an LP: band error of the returned knot vector g (both one-sided values of Γ_n at every knot) + repair cost + rounding allowance + 1e-9.
  - C_[a,b] = continuous functions on [0, 1], convex on [0, b], concave on [a, 1]. [a, b] lies inside one spacing, so no knot is strictly inside (a, b).
  - Knots: {0, every sample point x_(i), a, b, 1} (duplicates removed), so every jump of Γ_n is a knot. On each open knot interval g is affine and Γ_n is constant, so the supremum of |g − Γ_n| there is determined by the endpoint limits, which are the two one-sided values checked at the knots.
  - Let J_k = slope change of g at interior knot index k (knot indices 0 … K−1; ia, ib = indices of a, b). Convex side: k ≤ ib − 1 (so t_k ≤ a; LP rows j + 2 ≤ ib with k = j + 1). There δ_k = (−J_k)⁺, and we add δ_k(t_k − s)⁺, which is convex and zero on [a, 1]. Concave side: k ≥ ia + 1 (so t_k ≥ b; LP rows j ≥ ia). There δ′_k = (J_k)⁺, and we add −δ′_k(s − t_k)⁺, which is concave and zero on [0, b]. If a = b (ia = ib), the seam knot itself carries no constraint (no cross-seam slope condition), as in the class definition. If a < b, then ib = ia + 1 and the knot a is convex-side and the knot b is concave-side.
  - So the repaired g + h is in C_[a,b], and ‖h‖ ≤ Σδ + Σδ′ (the repair cost).
  - Rounding allowance A = 16u(Σ|slopes| + K·repair + 1) with u = 2⁻⁵³, plus the separate term 1e-9 (Lemma R below).
  - A non-finite value counts as an LP failure.
- LP failure after all fallbacks: inner node → bound +∞ (split further); leaf → Z = 1 at every cutoff (listed).
- Nested cutoffs: a higher cutoff is evaluated only after Z = 1 at the lower one.
- Not interval-checked: the uniform generator (numpy PCG64 via `default_rng(SeedSequence(...))`), HiGHS (scipy 1.17.1), numpy 2.3.5. Stated, not hidden.
- Replay: every zero classification can be reproduced deterministically from (root, n, shard, chunk, index) with the frozen code.

**Lemma R (floating-point allowance).** Model: IEEE-754 binary64, round to nearest, fl(x ∘ y) = (x ∘ y)(1 + δ) with |δ| ≤ u for + − × ÷ and √ (no overflow: |slopes| ≤ 1/min spacing < 1e15, asserted by finiteness; underflow only adds absolute errors ≤ 2⁻¹⁰⁷⁴ per operation, which the term 16u · 1 dominates). The returned g_k and t_k are exact dyadic numbers, and the bound is about the exact piecewise-linear function through (t_k, g_k).
- Band: each computed difference hi_k − g_k, g_k − lo_k has relative error ≤ u, and all values are in [−2, 2]. So the true band ≤ computed band + 2u.
- Slopes: m̃ = fl(fl(g_{k+1} − g_k)/fl(t_{k+1} − t_k)) = m(1 + θ), |θ| ≤ γ₃ = 3u/(1 − 3u). Jumps: J̃ = fl(m̃_{k+1} − m̃_k), so |J − J̃| ≤ γ₃(|m_k| + |m_{k+1}|) + u|J̃| ≤ 4.01u(|m̃_k| + |m̃_{k+1}|).
- (x)⁺ is 1-Lipschitz, so the true repair cost R ≤ Σ(−J̃)⁺ + Σ(J̃)⁺ + Σ|J − J̃| ≤ fl-sum·(1 + γ_K) + 8.02u·Σ|m̃|, because each slope enters two jumps.
- Hence true distance ≤ computed band + computed repair + 8.02uΣ|m̃| + 1.01Ku·repair + 2u ≤ band + repair + A − (a part of the "+1" term). The final float additions add ≤ 4u·ub < 1e-15.
- Threshold: t̃ = fl(c/fl(√n)) = (c/√n)(1 + θ₂). If the computed ub = band + repair + A + 1e-9 ≤ t̃, then the true distance ≤ t̃ − 1e-9 + 1e-15 ≤ c/√n · (1 + 2.01u) − 0.99e-9 < c/√n, since c/√n < 1. So every Z = 0 is backed by a true distance ≤ c/√n, i.e. Y = 0. ∎

## Part A — new certificate (fresh root)
- **Claim tested.** For n ∈ {800, 1600}: max_j P(√n S̄_{n,I_j} > c) ≤ α. The blocks I_j are as in DT8: [0.5, 0.9], [0.9, 0.97], [0.97, 0.995], [0.995, 1]. By Theorem U1 (DT7-UNIF, PROVEN), a pass gives sup over all continuous unimodal F of P_F(√n D_n > c) ≤ α.
  - α = 0.05: the single cutoff 0.7315.
  - α = 0.01: fixed sequence 0.87, then 0.86. At a look, 0.86 is certified only if 0.87 also passes at that look.
- **Root 2026110722** (fresh; CONFIRM). Shards s = 0…10. Chunk k of shard s uses `SeedSequence([2026110722, n, s, k])` with 500 samples (`np.sort(rng.uniform(size=n))` each); one JSON per chunk.
- **Looks** (chunk prefix 0…L−1 of EVERY shard; never "whichever finished first"):
  - n = 800: L = 36, 72, 108 → N = 198,000; 396,000; 594,000.
  - n = 1600: L = 24, 48, 72 → N = 132,000; 264,000; 396,000.
- **Decision at a look.** b = max{x : P(Bin(N, α) ≤ x) ≤ D} (−1 if none; exact recurrence, `accept_count` of dt8_cert_analysis.py). A cutoff passes iff all 4 block counts ≤ b.
- **Error budget.** D = 0.01/12 (2 n × 2 α × 3 looks).
  - Intersection–union: under H0 (max_j p_j ≥ α), pick a block j\* with p_{j\*} ≥ α. A false all-blocks certificate at a look needs block j\* to pass there, which has probability ≤ D.
  - The 1% fixed sequence adds no factor (the first false cutoff must pass).
  - So P_MC(any Part A claim false) ≤ 12·D = 0.01. The stopping look is one of the 3 covered looks, so reporting at it adds nothing.
- **Stop and report rules.**
  - An n stops after the first look at which 0.7315 (α = 0.05) is certified.
  - The 1% result is reported by the fixed sequence at that same look; no other look is substituted.
  - Not certified at look 3 → "not certified (inconclusive)", never "invalid".
  - Point estimates are reported, but they carry no certificate.
- **Power** (edge block alone, exact binomial; planning only):
  - n = 800, if the true rate is 4.86%: 0.39 / 0.82 / 0.97 by looks 1 / 2 / 3.
  - n = 1600, if 4.73%: about 0.8 / 0.99 / 1.0; lower if the true rate is nearer 4.9%.
- **Commands.** `python code/dt9_orchestrate.py 2026110722`, then `python code/dt9_cert_analysis.py 2026110722`. On a crash: rerun the same command (existing chunks are skipped; same seeds).
- **Time.** ≤ 3.2 h (n = 800) + ≤ 5.7 h (n = 1600) on 11 cores.

## Part R — certified replay of DT8 Part B (no new draws)
- Regenerate exactly the DT8 samples (root 2026110711; `SeedSequence([2026110711, n, shard])`; DT8 draw order) for n = 50, 100, 200, 400, 800, 1600 (N = 100k, 100k, 100k, 60k, 30k, 15k).
- **Identity check:** DT8 stored the edge value V for EVERY sample. Every replayed V must equal the stored V of the same sample within 1e-9, and the per-shard counts must match. Same generator and order as DT8 (`np.random.default_rng(np.random.SeedSequence([2026110711, n, shard]))`; `for _ in range(shard, R, 11): x = np.sort(rng.uniform(size=n))`; numpy 2.3.5).
- Compute the conservative indicators (above) at every DT8 grid cutoff (0.7315, 0.735, 0.74, 0.75, 0.86, 0.87, 0.88, 0.90).
- Apply the **unchanged DT8 rule**: d = 0.01/48; test from the largest cutoff down; stop at the first failure; all blocks must pass.
- A DT8 claim is **KEPT** iff the replay certifies a cutoff ≤ the DT8 one. Otherwise it is **WITHDRAWN**, and the replay result replaces it. The replay can never strengthen a DT8 claim.
- **Commands.** `python code/dt9_dt8_replay.py n R s 11` for every n and s = 0…10, then `python code/dt9_replay_analysis.py`. Run after Part A, or on spare cores; it uses no random draws, so its timing does not matter for validity.
- **Family error.** If every DT8 claim is KEPT, then the DT8 claims (≤ 0.01, from the DT8 rule and the replayed indicators) and the Part A claims (≤ 0.01) together have error ≤ 0.02. Before Part R finishes, no combined statement is made.

## Not claimed
Values of n outside the tables; that a non-certified cutoff is invalid; optimality of a cutoff; tied or rounded data; anything asymptotic (DT10).

## Code freeze (sha256, first 8; read 2026-11-07 01:56:16Z)
- code/dt9_block_cert.py cd8eaa67; code/dt9_cert_analysis.py 397b8d53; code/dt9_orchestrate.py 6fe0d681; code/dt9_dt8_replay.py 3f90085a; code/dt9_replay_analysis.py a68694c1.
- Imported and unchanged: code/sup_pinned.py 2f2441b1 (`_knots`, `pieces`, `SAFETY`, `TRIES`); code/pinned_dip.py 791dff28 (`V_left`); code/dt8_cert_analysis.py bf02a7ff (`accept_count`); DT7-UNIF.md f5f8825d.
- No change to any of these files after freezing. A crash is resumed with the same command and root; never fresh seeds.
