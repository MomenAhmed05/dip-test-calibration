# DT1b-PREREG — Step 3b confirmation. CONFIRM tier.

Frozen 2026-11-05 (file last modified 19:02:38 UTC), before any 3b output (result file last modified 19:15:06 UTC). Conventions: CONVENTIONS.md (CV-1..CV-10; roots in CV-6).

## 1. Statistic and calibration
- D_alg (CV-1: diptest 0.12, `allow_zero=False`; D := 0 for an all-identical sample). Reject iff D ≥ c (CV-4; no reference ties at q̂ in any earlier run; ties are counted again here).
- Uniform reference: R = 1,000,000 per n, root **2026110601**, 20 chunks of 50,000 (spawn (0, n, chunk)).
- n grid: {10, 20, 30, 50, 100, 200, 500, 1000, 2000, 5000, 12800}.
- Bracket: order statistics X_(l) ≤ q ≤ X_(u) with binomial tail ≤ γ/2 each. γ = 0.001 is split equally over the 22 (n, α) thresholds.

## 2. Laws and cells (100,000 datasets per cell; root **2026110602**, 10 chunks of 10,000, spawn (1, f, n, chunk))
| f | law | n grid | role |
|---|---|---|---|
| 0 | U(0,1) | all 11 n | anchor |
| 1 | W(0.1, 0.005) | 50 ... 12800 (8 n) | **primary** |
| 2 | ramp (w = 0.1, ε_r = 0.01; CV-9) | 50 ... 12800 | **primary** |
| 3 | W(0.4, 0.005) | 50 ... 12800 | secondary |
| 4 | probe step: density ∝ 0.033 on (0, 0.994), 1 on (0.994, 1) (the probe's printed, rounded "step#0") | 10 ... 200 (6 n) | probe re-derivation |
| 5 | probe edge spike: 0.9·U(0,1) + 0.1·(1 − \|N(0, 0.003²)\| mod 1) | 10 ... 200 | probe re-derivation |
| 6 | probe atom: 0.9·U(0,1) + 0.1·δ_1 | 10 ... 200 | probe re-derivation (atom row) |

All are IN (CV-10: R8 for f1-f4; f5 is a uniform plus a density decreasing away from 1; f6 is R3-type with an atom at the mode).

## 3. Primary endpoint (frozen statement)
**Claim E-3b:** for each of the 12 cells {f1, f2} × n ∈ {200, 500, 1000, 2000, 5000, 12800}, the rejection probability of the uniform-calibrated 5% dip test, at the own-n uniform quantile, exceeds 0.15.
- **Test:** an intersection-union test. Each cell is a one-sided exact binomial test at level δ = 0.001 against r ≤ 0.15, using the count of D ≥ X_(u), the upper bracket end. Equivalent: the one-sided Clopper-Pearson lower bound at 1 − δ is > 0.15.
- The claim holds iff all 12 cells pass. The overall error probability is ≤ δ + γ = 0.002.
- **Why 15%:** three times the nominal level. A test with this size behaves like a 15% test while it reports 5%.
- This is a claim about these 12 cells, not about "every n ≥ 200".

## 4. Secondary (reported, not part of the claim)
- **S1.** All cells at α = 5% and 1%: point rejection, plus simultaneous Bonferroni bounds (two-sided, level 1 − 0.05/K with K the number of reported cells) at X_(u) / X_(l).
- **S2.** Uniform anchor: the CP interval at q̂ contains α at every n (Bonferroni over 22 checks).
- **S3. Probe re-derivation table:** the probe value, the new value, the difference, and an explanation (A-2).
  - Probe values (`prior/exploratory/dip-probe/RESULTS.md`):
    - step 5%: 9.1 / 11.4 / 14.8 / 19.2 / 24.3 / 24.3;
    - step 1%: 1.9 / 3.5 / 4.4 / 6.5 / 9.8 / 10.5;
    - edge spike 5%: 7.4 / 8.7 / 10.8 / 14.2 / 21.2 / 26.1 (n = 10 / 20 / 30 / 50 / 100 / 200);
    - atom at n = 200: 26.5% (5%), 11.1% (1%).
  - **Agreement rule:** |new − probe| ≤ 1.0 point, or explained.
- **S4. Limit implications at n = 12800 (hybrid):** for f1, f2, f3, the five-point check of Z_n = 2√n·D/√(1−w) against F_M (f1, f2) and against F_M·H_U (f3, with H_U from this run's n = 12800 reference). Simultaneous 95% CIs (Bonferroni over 15) must lie inside ±0.01 for "five-point equivalence"; otherwise "not shown". Also the predicted rejection 1 − F_M(2q̂/√(1−w))·H_U(q̂/√w), using this run's q̂ at n = 12800, against the observed rejection.
- **S5.** S_n = √n·[D − max(V_sh, T_top)] and √n·(D − V_sh) at n = 12800 for f1, f2, f3: quantiles (1%, 50%, 99%, max), and P(D < V_sh − 1e-12).
- **S6. LP audit (independent implementation):**
  - random: 100 samples per (f, n) for f ∈ {0, 1, 2, 3} and n ∈ {50, 100, 200}, chosen by root 2026110603;
  - targeted: every sample of those cells with |D − q̂_.05| < 1e-6.
  - Rule: |D_alg − D_LP| ≤ 1e-9 for all audited samples. Any failure → audit and report before reading S1-S5.

## 5. What would falsify
- Any primary cell fails → claim E-3b is not confirmed (report the failing cells).
- The LP audit fails → stop; audit the implementation.
- The uniform anchor fails S2 → the reference is suspect; the claim is not read.
- Probe re-derivation off by > 1 point without an explanation → report as an unexplained mismatch.

## 6. Recorded S predictions (before running)
- E-3b: all 12 cells pass. W(0.1) is about 25-26% at n = 200-2000 and about 24-25% at 5000-12800; the ramp is about 25%.
- W(0.4): about 7-8% at n ≥ 1000.
- Probe re-derivation: all within 1 point (the step law's heights are rounded in the printout, so it may differ by up to ~1 point).
- S4:
  - W(0.1) and the ramp: five-point equivalence shown at 12800 (the 3a max error was ≤ 0.0092);
  - predicted vs observed rejection within 0.5 point.
- S5: P(D < V_sh) = 0 for f1-f3.
- S6: all audits pass.
