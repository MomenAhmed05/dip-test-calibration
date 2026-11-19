# DT2-PREREG — Step 4 conjecture farm, part 1 (block formula). DEV tier.

Frozen 2026-11-05 (file last modified 19:25:26 UTC), before any DT2 output (result file last modified 19:34:42 UTC). Root **2026110610** (DEV; declared in CV-6).

## 1. Conjecture under test (BF = block formula)
For a unimodal step density with maximal flat blocks i (masses p_i), and a modal block m (the unique highest block, mass p_m):

√n·D_n ⇒ L = max( max_{i ≠ m} √p_i·M_i/2 , √p_m·U′ ),

with all terms independent; M ~ Balabdaoui-Pitman; U′ ~ the uniform dip limit. The predicted rejection at level α is

r_BF = 1 − Π_{i≠m} F_M(2q/√p_i) · H_U(q/√p_m),

where q = √n·q̂_n(α) and H_U = the CDF of √n·D_n under the uniform, both from **this run's own-n reference** (the hybrid finite-n version). F_M is from the B-F 2011 Table 2 (x ∈ [0.90, 1.79], linear interpolation; F_M := 1 above 3.0 and a table lookup otherwise; if 2q/√p_i > 1.79, the B-F Table 5 tail value is used, F_M(1.79) = 0.99386, and F_M(x) := 1 − 0.00614·exp(−8(x − 1.79)) beyond; recorded as approximate).

## 2. Laws (density heights h on intervals; all unimodal by construction: nondecreasing to the modal block, then nonincreasing)
| f | name | blocks: interval → mass | role |
|---|---|---|---|
| 0 | W(0.1, 0.005) | (0, .995) → .9; (.995, 1) → .1 [modal] | anchor |
| 1 | two shoulders, big jump | (0, .5) → .3; (.5, .995) → .6; (.995, 1) → .1 [m] | join |
| 2 | two shoulders, small jump (heights 0.86, 0.9495) | (0, .5) → .43; (.5, .995) → .47; (.995, 1) → .1 [m] | join / slow merger |
| 3 | three-block staircase | (0, .33) → .2; (.33, .66) → .3; (.66, .995) → .4; (.995, 1) → .1 [m] | join |
| 4 | interior plateau, unequal | (0, .6) → .6; (.6, .605) → .1 [m]; (.605, 1) → .3 | modal + two sides |
| 5 | interior plateau, symmetric | (0, .4975) → .45; (.4975, .5025) → .1 [m]; (.5025, 1) → .45 | modal + two sides |
| 6 | W, narrow top | (0, .9995) → .9; (.9995, 1) → .1 [m] | geometry (same BF limit as f0) |
| 7 | W, wide top | (0, .95) → .9; (.95, 1) → .1 [m] | geometry (same BF limit as f0) |
| 8 | W reflected | (0, .005) → .1 [m]; (.005, 1) → .9 | reflection |
| 9 | shrinking jump | as f2 but heights h1 = c, h2 = c(1 + 4/√n), masses solved so that the total on (0, .995) is .9 | triangular array |
| 10-13 | K4 path ε = 0.25w | w ∈ {.05, .02, .01, .005}: (0, 1 − ε) → 1 − w; (1 − ε, 1) → w [m] | sup approach |
| 14 | W(0.02, 0.005) | (0, .995) → .98; (.995, 1) → .02 [m] | the 9.7% puzzle trajectory |

n:
- f0-f9: {800, 3200, 12800};
- f10-f14: {3200, 12800, 51200}.

20,000 datasets per cell.

Reference: R = 1e6 for n ≤ 12800 and 5e5 at n = 51200. Spawn (0, n, chunk); laws (1, f, n, chunk).

## 3. Per-sample quantities
- D (D_alg).
- The block terms:
  - rising-side block i: V_i = the convex distance of F_n on the block's points (`vn_sample` with its rank offset);
  - falling-side block: the concave distance (the same function on the reflected points);
  - modal block: T_m = (k_m/n)·D_alg(block points).
- BF_n = max of the block terms. R_n = √n(D − BF_n). Report q01 / q50 / q99 / max, and P(|R_n| ≤ 1e-12).

## 4. Readings fixed in advance (α = 5%)
- **BF consistent at n** iff |obs − r_BF| ≤ 1.5 points (about 2.5 SE at 20k reps, plus finite-n allowance).
- **Primary farm reading P1:** BF is consistent at n = 12800 for every f ∈ {0, 1, 3, 4, 5, 6, 7, 8} → "BF SURVIVES (joins, modal, geometry, reflection)". Each failing f is reported as a kill of BF for that structure.
- **f2 (small fixed jump):** BF predicts the separated-block value r_BF(f2). The S prediction is a trajectory from near W(0.1)'s value (merged) at n = 800 towards r_BF(f2), not yet consistent at 12800. Reading: "slow merger" if obs decreases monotonically with n towards r_BF.
- **f9 (shrinking jump):** S predicts a roughly constant rejection in n, between r_BF(f2) and the merged value. Reading: "local transition" if the n = 3200 and 12800 values differ by ≤ 2 points.
- **K4 (f10-f13):** obs approaches r_BF(f) as n grows: |obs − r_BF| non-increasing in n and ≤ 2 points at 51200.
- **f14:** obs increases with n towards r_BF ≈ 32%: "pinning transition" if monotone increasing over 3200 → 51200.
- **R_n:** q99 of |R_n| decreasing in n for f0-f8.

## 5. Recorded S predictions
- P1: BF SURVIVES for all 8 structures. Largest misfit for f5 (interior symmetric: two shoulders + modal term all similar).
- r_BF values (approximate, n = 12800 q ≈ 0.543):
  - f0 ≈ 24.8%;
  - f1 ≈ 1 − F_M(1.40)·F_M(1.98) ≈ 7.5%;
  - f3: largest block .4 → ≈ 1 − F_M(1.72)·F_M(1.98)·F_M(2.43) ≈ 1.5%;
  - f4 ≈ 1 − F_M(1.40)·F_M(1.98) ≈ 7.5%;
  - f5 ≈ 3.5%;
  - f6 = f7 = f8 = f0 ≈ 24.8%.
- f2 slow merger (≈ 22% at 800 → ≈ 12% at 12800; r_BF ≈ 5%).
- K4: w = .05: ≈ 29% (obs a bit lower); w = .005 at 51200: about 25-30% (not yet at r_BF ≈ 31%).
- f14: 15% → 25% → 30% over 3200 → 51200.
