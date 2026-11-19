# DT8-PREREG — Step 10 Monte Carlo certificates (gaps G1 and G2)

Status: **FROZEN** 2026-11-06 22:24Z (`date -u`; code hashes read 22:23:26Z).

## Disclosure: what was seen before freezing (DEV only)
- Root 2026110701 (pinned Y(u)/Z(u), n = 200/1000/4000): tails increase in u; the edge u = 1 is the worst.
- Root 2026110703 (block pilot, `dt8_block_dev.py`, n = 50…800, 4k–20k reps): the worst block is always the edge block [0.995, 1]. Tails at c = 0.7315: 3.8%, 4.2%, 4.6%, 4.6% and 4.3% for n = 50, 100, 200, 400 and 800. At c = 0.86: 0.7–1.0%.
- Roots 2026110712/13: tiny quick runs of both drivers (code checks only).
- The pilot chose the blocks, the cutoff grids and the N values below. CONFIRM roots were never used before.

## Part B — finite-n uniform-validity table (gap G2; DT7-UNIF Theorem U1)
- **Claim tested.** At each listed n, the cutoff c satisfies max_j P(√n S̄_{n,I_j} > c) ≤ α. By Theorem U1 (PROVEN), the test "reject if √n D_n > c" then has size ≤ α **uniformly over all continuous unimodal F**.
- **Blocks** I_j: [0.5, 0.9], [0.9, 0.97], [0.97, 0.995], [0.995, 1]. Branch-and-bound floor 0.65 (√n scale); this is below every tested cutoff, so the exceedance indicators are exact.
- **n and N:** 50, 100, 200 → 100,000; 400 → 60,000; 800 → 30,000; 1600 → 15,000.
- **Cutoff grids**, tested from the largest down (fixed sequence, stop at the first failure):
  - α = 0.05: 0.75, 0.74, 0.735, 0.7315.
  - α = 0.01: 0.90, 0.88, 0.87, 0.86.
- **Decision.** A cutoff passes iff, for all 4 blocks, the count K_j ≤ b(N, α, d) = max{x : P(Bin(N, α) ≤ x) ≤ d}, with d = 0.01/48. This is the same as the one-sided Clopper–Pearson upper limit at error d being ≤ α. The counts are computed with exact rational-recurrence arithmetic (mpmath, 50 digits). Reported: the smallest cutoff that passed before the first failure, or "none".
- **Failure policy.** A block whose LP fails after all fallbacks counts as an exceedance (+∞). It is not dropped.
- **Guarantee.** P_MC(any reported cutoff is false) ≤ 48·d = 0.01 (Bonferroni over 12 (n, α) pairs × 4 blocks); a fixed-sequence argument gives the sharper bound 12·d = 0.0025. Samples are shared across blocks, cutoffs and α; the fixed-sequence argument covers this.
- **Root** 2026110711, shards s = 0…10 of 11 (`SeedSequence([root, n, shard])`). Command: `python code/dt8_block_dev.py n N s 11 2026110711 cert`.
- **Not claimed.** Values of n outside the table. That a non-certified smaller cutoff is invalid. That the certified cutoff is optimal.

## Part A — certified bracket for q̄_α (gap G1 remainder)
- **Object.** q̄_α = (1 − α) quantile of D(B), the H&H dip functional of a standard Brownian bridge (the limit of √n·D_n under the uniform; H&H Thm 3).
- **Per path** (`code/qbar_cert.py`):
  - Exact Gaussian skeleton on m = 2,000,000 cells.
  - D(B_m) = λ·dip(ECDF of N = 8,000,000 midpoint quantiles of H = id + B_m/λ), λ = 1.05·max|slope|. Justification: H&H Thm 1 and Thm 2, affine and scale invariance; this is not a cone claim.
  - E_m drawn from P(E_m ≤ x) = K(x√m)^m (K = law of sup|bridge|).
  - lo = max(0, λ(dip − η) − E_m) and hi = λ(dip + η) + E_m, with η = 1/(2N) + 1e-12.
  - Coupling statement: the skeleton and E_m admit a coupling with a standard bridge B under which lo ≤ D(B) ≤ hi almost surely.
- **R = 12,000 paths.** Root 2026110710, 11 shards. Command: `python code/qbar_cert.py 2000000 8000000 12000 s 11 2026110710 acert`.
- **Bracket.** Upper bound hi_(k), with k the smallest index such that P(Bin(R, 1 − α) ≥ k) ≤ ε. Lower bound lo_(ℓ), with ℓ the largest index such that P(Bin(R, 1 − α) ≤ ℓ − 1) ≤ ε. Here ε = 0.01/4 and α ∈ {0.05, 0.01}.
- **Worst-case size** of the uniform-calibrated test, 1 − F_M(2q̄_α) (DT3 Thm 1′/Thm 4, C4). It is bracketed with the certified F_M enclosure `deposits/dt8/fm_cert_h0.0002.npz`, using outward grid points.
- **Guarantee.** P_MC(both brackets hold for both α) ≥ 0.99. Combined with Part B: ≥ 0.98 (union bound).

## Numerical caveats (stated, not hidden)
- **Floating point is not interval-checked.** This covers the Gaussian generator, the LP solver, diptest 0.12.0 (Maechler's corrected AS 217), the Kolmogorov inverse (scipy 1.17.1) and numpy 2.3.5.
- **Safety margins.** LP values get +1e-9 and the realised band error; the worst slope violation is logged (smoke test: 9e-16). ECDF slack is +1e-12. The rigour claim is: "Monte Carlo certificate with simultaneous coverage ≥ 0.98, under the ideal arithmetic and iid-uniform model, with the stated floating-point margins." It is not a deterministic proof.
- **diptest exactness check.** It agrees with our independent LP dip to 3.6e-13 on 20 samples (`code/pinned_dip.py` self-test).

## Code freeze (sha256, first 8; read 2026-11-06 22:23:26Z)
- code/sup_pinned.py 2f2441b1; code/dt8_block_dev.py de7f9c27; code/qbar_cert.py 9506890b; code/dt8_cert_analysis.py bf02a7ff; code/pinned_dip.py 791dff28; DT7-UNIF.md a2a96069; deposits/dt8/fm_cert_h0.0002.npz 943658cc.
- Analysis commands after the runs: `python code/dt8_cert_analysis.py B 2026110711` and `python code/dt8_cert_analysis.py A 2026110710 acert`.
- No change to any of these files after freezing. If a run crashes, resume the same shard with the same root; never draw fresh seeds.
