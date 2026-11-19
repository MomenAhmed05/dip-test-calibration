# DT0d-PREREG — Step 2c follow-up screen. DEV tier only.

Frozen 2026-11-05 (file last modified 18:34:46 UTC), before any DT0-d output existed (result file last modified 18:37:51 UTC). Origin: `DT0c-SCREEN.md` §6.

**Tier:** DEV (exploratory, never confirmatory). Nothing here changes the frozen DT0-c verdict.

## 1. Statistic, calibration, seeds
- D_alg (CV-1). Own uniform reference H_n, R = 1,000,000 per n, for n ∈ {200, 800, 3200, 12800}.
- Reject iff D ≥ c (CV-4). α = 0.05 is primary; α = 0.01 is reported.
- Bracket γ = 0.001 over the 8 (n, α) pairs; bounds at X_(u) / X_(l) as in DT0-c.
- Seeds: new DEV root **2026110550** (declared in CV-6 before the run). Reference spawn (0, n, chunk); laws spawn (1, f, n).
- 20,000 datasets per cell. Every dip is saved to `deposits/dt0d/`.

## 2. Laws (f index; membership certificates per CV-10)
**(i) Local slope scaling.** Witness W(0.1, 0.005) with shoulder density (0.9/L)[1 + β(x/L − ½)]:
- f = 0..3: β_n = c/√n with c ∈ {0, 2, 4, 8}. c = 0 is W itself (the anchor).
- f = 4, 5: fixed β = 0.1 and 0.25 (comparators).
- IN by R8 (the density is nondecreasing; the shoulder end density is ≤ 1.81 < 20).
- Note: β_n ≤ 2 for all cells (max is c = 8 at n = 200: β = 0.566).

**(ii) Fixed Gaussian noise on W.** f = 6, 7, 8: W + σE with σ ∈ {0.005, 0.025, 0.05}. IN by R7.

**(iii) Laws with a flat shoulder for a structural reason:**
- f = 9..12 — **p-value mixture** π0·U(0,1) + (1 − π0)·Beta(a, 1), with (π0, a) ∈ {(0.9, 0.1), (0.9, 0.25), (0.8, 0.1), (0.8, 0.25)}. Density π0 + (1 − π0)·a·x^{a−1} on (0, 1) is strictly decreasing, so the law is IN with mode 0 (R4 type, bounded support).
- f = 13, 14 — **one-sided z-test p-values** π0·U + (1 − π0)·law of 1 − Φ(Z + μ), with π0 = 0.9 and μ ∈ {2, 3}.
  - Alternative density φ(Φ⁻¹(1−p) − μ)/φ(Φ⁻¹(1−p)) = exp(μ·Φ⁻¹(1−p) − μ²/2). This is decreasing in p for μ > 0.
  - So the mixture density is decreasing: IN, with mode 0.
- f = 15, 16 — **biased-forecast PIT** u = Φ(Z − δ), δ ∈ {0.25, 0.5}.
  - Density exp(−δ·Φ⁻¹(u) − δ²/2) is decreasing in u: IN, with mode 0.
  - This is a single-edge pile-up. Two-edge PIT (dispersion errors) is NOT unimodal and is not used.
- f = 17 — anchor **U(0,1)**.

## 3. Reported quantities
For every cell: the point rejection at q̂ (α = 5%, 1%); 95% CP bounds at X_(u) / X_(l); and the quantiles (0.5, 0.9, 0.95) of √n·D_n.

## 4. Readings fixed in advance (descriptive rules, α = 5%)
- **Anchors:** U must be CP-consistent with α at all n (Bonferroni over 8 checks, as in DT0-c). W (c = 0) must be ≥ 15% at n = 200. Otherwise the screen is VOID.
- **(i) Local scale supported** iff all three hold:
  - (a) for each c, the rejection at n = 3200 and n = 12800 differ by ≤ 2.5 points;
  - (b) at n = 12800 the rejection is nonincreasing in c;
  - (c) fixed β = 0.1 at n = 12800 is below its own value at n = 800 by > 2.5 points.
  - Otherwise: NOT SUPPORTED, with the pattern reported.
- **(ii) Turnover seen** for a σ iff the 95% CP upper bound at n = 12800 is below the 95% CP lower bound at n = 3200. Otherwise: no turnover up to n = 12800.
- **(iii) MATERIAL** iff the 95% CP lower bound (one-sided, Bonferroni over the 8 laws × 4 n = 32 cells, i.e. level 1 − 0.05/32) is ≥ 7.5% at some n. CONTROLLED iff every upper bound is < 7.5%. Otherwise INCONCLUSIVE.

## 5. Recorded S predictions (before running)
- (i):
  - c = 0 stays at about 25%;
  - c = 2 is about 15-20%, roughly stable in n;
  - c = 8 is ≤ 7%;
  - fixed β = 0.1 falls below 7.5% by n = 12800;
  - local scale SUPPORTED.
- (ii): no turnover up to n = 12800 for σ = 0.005 and 0.025. For σ = 0.05: no turnover, but a lower level (the spike is smeared over about 0.05).
- (iii):
  - p-value mixtures with a = 0.1: MATERIAL (≥ 15% at n ≤ 800), decaying slowly at large n;
  - a = 0.25: smaller, possibly INCONCLUSIVE;
  - z-test p-values: μ = 3 MATERIAL, μ = 2 weaker;
  - PIT tilt: CONTROLLED (the whole density is tilted, so no flat stretch).

## 6. What this screen does NOT establish
- A size.
- A limit.
- Real-data prevalence: these are named laws.
- A real use of the dip on p-values or PIT values (H-002 must find one).
