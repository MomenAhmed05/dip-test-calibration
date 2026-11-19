# DT3c-PREREG — Step 5 CONFIRM: the block formula (BF) at large n + pinning-overshoot check

Frozen 2026-11-05 23:36Z (`date -u`). Quick dry run on DEV roots 2026110632/33 only (9 s, code paths checked; CONFIRM roots unseen). Tier: **CONFIRM** (fresh roots, frozen before any run).

## Purpose
1. **BF (Conjecture C1, Theorem 2(b), Lemma G).** Test √n·D_n ⇒ max_i(√p_i·M_i/2 ; √p_m·U′) on five structurally different FIXED step laws, at n up to 204,800, on the whole distribution (five tail points), and on the paired remainder R_n = √n·(D_n − BF_n), where BF_n = max of the finite-n block functionals (non-modal blocks: V on the block; modal block: (k/n)·dip(block sub-sample)).
2. **Pinning overshoot (C3/C4 benchmark).** Re-test the four hardest DEV pinning cells at larger n against the benchmark b(n) = 1 − F_M(2·q_n), q_n = own-n uniform reference quantile. An overshoot refutes the BENCHMARK as an asymptotic upper bound; it does not by itself refute R_unif = R_pt.

## Objects (CONVENTIONS v1)
- D_alg = diptest dip, allow_zero = False (`dt0c_screen.dips_of`); V via `path_smoke.vn_sample`; block terms via `dt2_farm.block_terms`; F_M via `dt1b_confirm.bf_table` (Balabdaoui-Filali grid + interpolation; numerical, Step 6 certifies).
- Rejection: √n·D_n > q_n(α) (strict), q_n = order-statistic quantile of the own-n uniform reference (method "inverted_cdf").

## Laws (fixed in n; densities shown)
| id | edges | masses | modal block | densities | BF 5% / 1% (pred., DEV ref n=12,800) |
|---|---|---|---|---|---|
| L1 interior modal plateau, blocks on both sides | 0, .30, .45, .55, 1 | .24, .15, .20, .41 | 2 | .80, 1.0, 2.0, .911 | 1.19 / 0.15 |
| L2 boundary modal block = W(0.1, .005) (anchor) | 0, .995, 1 | .90, .10 | 1 | .905, 20 | 24.79 / 10.49 |
| L3 strongly unbalanced masses | 0, .98, .99, 1 | .95, .02, .03 | 2 | .969, 2, 3 | 28.22 / 12.55 |
| L4 close-but-unequal adjacent heights | 0, .50, .995, 1 | .43, .47, .10 | 2 | .86, .949, 20 | 3.91 / 0.64 |
| L5 several joins (5-step staircase) | 0, .2, .4, .6, .8, 1 | .12, .16, .20, .24, .28 | 4 | .6, .8, 1.0, 1.2, 1.4 | 0.02 / 0.00 |

Pinning cells (DT2b `spec`: top width ε = λ/√n, excess MASS e = κ/√n, w = ε + e(1 − ε)): P1 (κ 1.4, λ .1), P2 (1.7, .1), P3 (1.2, .1), P4 (2.5, 1). DEV at n = 12,800: 30.64 / 30.51 / 30.03 / 29.33% vs b = 31.70%.

## Design
- n ∈ {3,200; 12,800; 51,200; 204,800}. **Nested common random numbers:** each replicate draws one stream of 204,800 uniforms u; for each n, the sample is the prefix u[:n] mapped through that law's quantile map at that n (`dt2_farm.sample` logic). So all n share randomness (CRN across n and, for the reference, nested prefixes).
- Reference: 1,000,000 uniform replicates (nested prefixes give all four n). Root **2026110630**, spawn (0, chunk).
- Laws: 100,000 replicates per law (all four n from the same replicate). Root **2026110631**, spawn (1, law_index, chunk). law_index 0-4 = L1-L5, 5-8 = P1-P4.
- MC SE at 25-30% ≈ 0.14-0.145 pt (independent of the 1M-reference error).

## Readouts (all deposited in `deposits/dt3c/confirm.json` + per-cell .npy)
For each law and n: obs rejection at 5%/1%; r_BF(n) (formula with own-n q and own-n reference H for the modal block); five tail points t ∈ ref quantiles {.50, .75, .90, .95, .99}: P̂(√n D > t) vs BF prediction 1 − Π_i[...]; R_n quantiles {.5, .99, .999, max}, P(|R_n| > 0.01), P(|R_n| > 0.001), P(R_n < −1e-9) (D below BF_n would contradict the lower bound); for pinning cells: obs vs b(n), with SE.

## Decision rules (frozen)
- **BF CONFIRMED for a law** if at n = 204,800: |obs − r_BF| ≤ 1.0 pt (5%) and ≤ 0.7 pt (1%), all five tail-point differences ≤ 1.5 pt, AND q99|R_n| at 204,800 ≤ 0.5 × its value at 3,200 (or ≤ 0.002).
- **BF REFUTED for a law** if at both 51,200 and 204,800 the 5% (or any tail-point) difference exceeds 1.5 pt with the same sign and does not shrink (|diff₂₀₄₈₀₀| ≥ |diff₅₁₂₀₀| − 0.3).
- Otherwise **INCONCLUSIVE (slow)**.
- **Lower-bound sanity:** P(R_n < −1e-9) must be 0 for L2 and L4 when the shoulder has ≥ 2 points (Theorem 2(a) is a lower bound w.p. → 1; isolated small-n violations are recorded, not decisive).
- **Pinning:** OVERSHOOT if obs > b(n) + 3·SE at both 51,200 and 204,800; NO-OVERSHOOT if obs ≤ b(n) + 1·SE at all n; else UNCLEAR.

## Predictions (S, before running)
- L1-L5 BF CONFIRMED; L4 slowest (DEV f2: excess above BF shrinking with n); L5's 5%/1% rates are ≈ 0 so its test rests on tail points and R_n.
- R_n: P(|R_n| > 0.01) falls with n for every law (Lemma G).
- Pinning: NO-OVERSHOOT; obs stays 29-31% at all n.

## Not claimed by this run
No theorem is proved by simulation. A CONFIRMED BF supports Conjecture C1 / Lemma G only. F_M values are numerical (Step 6).

## Compute
Estimated 1.5-2 h on 11 workers (≈ 30 ms per replicate at n = 204,800). This run decides between "right limit with finite-n error" and "slightly wrong limit", which a single n cannot.
