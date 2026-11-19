# DT5b-PREREG — post-DT5 SUPPLEMENT: Cheng & Hall (1998) calibrated dip test vs uniform and LF

Status: **FROZEN** 2026-11-06 ~16:12Z (`date -u`). Draft ~16:02Z, written AFTER DT5 (CONFIRM) and AFTER one DEV quick run (§5); revised after one review round: convolution-length bug fixed, uniform node d = 0 added, near-2π normal band, fidelity gates F1-F4 (run in DEV, §6), wording.
Code: `code/ch98.py` (calibration), `code/dt5b_ch98.py` (runner + summary). Laws and samplers from `code/dt5_power.py` (sha 1d0f2c52); dip from `code/dt0c_screen.py` (7fe20a4f).

## 1. The procedure under test (from the paper, Sections 2-3)
Normal-kernel f̂ and f̂'' with normal-reference AMISE bandwidths (h0 = 1.0592·S·n^(−1/5), h2 = 0.9397·S·n^(−1/9)), 512-point grid, x̂0 = argmax f̂, d̂ = |f̂''(x̂0)| / f̂(x̂0)³. Calibration law with the same d: Beta(β, β) if d < 2π, normal if d = 2π, Student t_ν if d > 2π. Reject if dip > (1 − α) quantile of the dip under that law at the same n.
Our implementation choice (stated): the per-dataset bootstrap is replaced by a pre-computed table of dip quantiles at node d = 0 (uniform = Beta(1,1)) and 49 log-spaced d values in [10⁻³, 10³] (M = 20,000 samples per (n, node)); linear interpolation in d on [0, 10⁻³] and in log d above; clipping above 10³ (frequency reported); |d − 2π| < 0.1%·2π → normal. This APPROXIMATES the per-dataset bootstrap rule (Monte Carlo node error + interpolation error + clipping + KDE binning); it is not its exact large-B limit. The fidelity checks in §6 bound these errors.
Scope (paper's own assumptions): the theory assumes a unique mode with f''(x0) < 0; Ameijeiras-Alonso, Crujeiras & Rodríguez-Casal (2019, TEST 28:900-919) §2.3 state such regularity conditions. Laws with flat parts violate them. **We evaluate this implementation outside its regularity assumptions to determine its robustness on flat and pile-up distributions.** This is not a refutation of the regular-case theorem. Cheng & Hall (1999, AoS 27:1294-1315) already warn that the no-shoulder calibration tends to be anticonservative for a mode with a shoulder; our flats are a different object (positive-length flat parts), and atoms are not used here.

## 2. Design
- Same samples for uniform (DT5 own-n reference cutoff from `deposits/dt5/power.json`), LF (0.73139526 / 0.85987427 over √n) and CH98.
- Laws: the 9 DT5 null laws + 5 DT5 alternatives (sym 1.5, sym 2, asym(0.25, 2.5), valley(0.5, 1), satellite). n ∈ {50, 200, 800, 3200, 12800}; α ∈ {0.05, 0.01}. R = 10,000 per (law, n).
- Roots: table 2026110680, laws 2026110681 (fresh; banned after). DEV quick roots 2026110682/83 (seen).

## 3. Rules
- **S1 size:** per (method, law, n, α), two-sided CP-95%. Labels as DT5 R1 (tolerance 0.5 pt at 5%, 0.25 pt at 1%): certified within tolerance / demonstrated above tolerance / inconclusive, over all 9 null laws.
- **S2 implementation gates (all run in DEV before freezing, §6; re-reported from the CONFIRM run where applicable):** F1 binned vs unbinned d̂ relative difference ≤ 2%; F2 inversion round trip has 0 failures on [10⁻⁴, 10⁴] incl. the 2π band; F3 table interpolation vs direct quantile ≤ 1% relative at 6 between-node d; F4 reproduction of the paper's Table 1 (α = 0.05; normal, t6, beta(4,8), beta(12,4); n = 50, 200) within 3 pt of the paper's values (their SE ≈ 1 pt from 500 replications). Per-cell intervals are pointwise CP-95%; the S1 labels are per fixed method (intersection over cells), not simultaneous claims.
- **S3 power:** reported only for methods certified at that (n, α); otherwise the rejection rate is labelled "rejection rate (size not controlled)".

## 4. Predictions
- Q1 (DEV-informed): CH98 demonstrated above tolerance at every n, both α (uniform-law rate ≥ 25% at 5% for n ≥ 200).
- Q2 (DEV-informed): CH98 on N(0,1): within [2%, 8%] at 5% for all n.
- Q3: CH98 on W(0.1) and pin laws ≥ 75% at 5% for n ≥ 800.
- Q4: uniform and LF reproduce DT5 within CP intervals (sanity).

## 6. Fidelity gates (DEV, `code/ch98_fidelity.py`, root 2026110684, 16:07Z, `deposits/dt5b_fidelity.json`)
F1 max relative difference 0.83% (W(0.1), n = 800; all others ≤ 0.02%; two-point case fixed). F2 0 failures. F3 ≤ 0.7%. F4 at 5% (ours vs paper): normal 4.75 / 4.45% vs 3.6 / 4.8%; t6 2.95 / 2.70 vs 3.2 / 2.4; beta(4,8) 6.20 / 4.65 vs 5.6 / 4.2; beta(12,4) 5.10 / 5.30 vs 6.2 / 5.6 — all within 1.2 pt. Not in the paper: beta(2,2) 9.0 / 9.5% (a curved law; noted, not explained); uniform 29.8 / 48.9%. Extreme Student references (d = 50, 200, 1000; max |x| up to 7e45): diptest dip = LP dip = dip after an affine rescale (5/5 each).

## 5. DEV quick run seen — disclosure (before the bug fix; values change slightly after it)
`dt5b_ch98.py --quick` 15:52-15:54Z (`date -u`), roots 0682/83, M = 2,000, R = 1,000: CH98 at 5% — uniform law 28.6 / 46.1 / 61.8 / 74.1 / 82.2% (n = 50 … 12,800); W(0.1) 50-99.9%; pin laws 32-98.5%; normal 3.4-5.7%; beta(2,2) 7.1-10.8%. Alternatives: CH98 rejects more than uniform everywhere.
