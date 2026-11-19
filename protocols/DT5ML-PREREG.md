# DT5ML-PREREG — Step 7 ML items: O5-ML, gap statistic, O8 audit

Status: **FROZEN** 2026-11-06 ~15:45Z (`date -u`). Draft ~15:38Z, written AFTER DT5-RESULT (DT5 seen) and after the DEV runs listed in §6; revised after one review round (rule labels, K = 70, R = 50,000 for O8, 216-cell Bonferroni, baseline margin, adjudication, wording). Design source: DT5-ML-DESIGN v0.2.
Pipeline (hashes at freeze): `code/o5_ml.py`, `code/o5_gap.py`, `code/o8_audit.py`, frozen calibrator `deposits/o8/calibrator.pt` (trained 11:40-12:15Z, DEV root 2026110670, 1,500 iterations; `deposits/o8/train.json`), `code/fm_table.py` (25abe2ee), `code/dt0c_screen.py` (7fe20a4f). A separate Python environment (ClustPy 0.0.3) for o5_ml/o5_gap.

## 1. O5-ML (calibration-substitution sensitivity; DT5-ML-DESIGN §1)
- Methods and defaults: UniDip (1-D; α = 0.05 and 0.01) PRIMARY; SkinnyDip (d = 5; α = 0.05) extension; DipMeans (d = 5; n = 500; α = 0.001) secondary mechanism study. ClustPy defaults otherwise.
- Arms: uniform = ClustPy `dip_pval` (table); LF-asymptotic = 1 − F_M(2√n·D) with F_M from `fm_table.FM_pub`. Same data, same algorithm seeds.
- Laws (19): U, G, N1 = W(0.1), N2(σ) σ ∈ {0.0005, 0.002, 0.01, 0.05}, N2-pin(κ) κ ∈ {0.85, 1.4, 2.5}, R0(σ) σ ∈ {0.002, 0.05}, A_N2, A_R0 (atoms; reported separately), price laws ReLU(N(0.5,1)), mix(2.5/3/4), tanh(3·flat).
- n ∈ {500, 2000}; R = 2,000 per (method, law, n) for UniDip/SkinnyDip, 1,000 for DipMeans. Roots 2026110660 (data) / 2026110661 (algorithm init); spawn (root, law, n, rep).
- Outcome: P(k̂ > 1) per arm; Δ = P_U − P_LF; both discordant proportions; exact McNemar; Newcombe CI for Δ (paired).
- Rules (analyser `code/o5_analysis.py`, frozen): **M-R1** "positive calibration-substitution effect on over-splitting" is declared in a cell when the exact one-sided McNemar p = P(Bin(b + c, ½) ≥ b) < 0.05/K and b > c (b = uniform-only splits, c = LF-only; p = 1 when b + c = 0). Family (machine-enumerated by the analyser): the 10 non-control continuous null laws N1, N2(×4), N2-pin(×3), R0(×2) × {UniDip α 0.05/0.01 × 2 n, SkinnyDip 0.05 × 2 n, DipMeans n = 500} → **K = 70**. Δ CI: Newcombe (1998) method 10 (paired score interval). It identifies the effect of swapping the calibration inside the fixed pipeline; it does not attribute every removed split to the root test. **M-R2** U, G, atom laws and price laws are reported only. **M-R3** no pipeline-level size claim; "M-R1 not declared" is non-detection, not equivalence.

## 2. Gap statistic (secondary; not a calibration comparison)
ClustPy `GapStatistic` k = 1..6, B = 100, 1-SE rule; n ∈ {100, 500}; R = 200 (reduced from 1,000 at freeze for cost, stated); root 2026110665. Laws: N1, N2(0.002), R0(0.002), A_N2 (pile-ups) vs U, G, logistic (log-concave controls). Outcome P(k̂ > 1) with CP-95%. R = 200 outer datasets vs B = 100 inner reference sets are different randomness; k = 6 boundary handled as in ClustPy 0.0.3 (no change). The gap statistic targets a reference-based notion of clusters, not a level-α unimodality test: k̂ > 1 on a unimodal law is not a size violation (Tibshirani, Walther & Hastie 2001). Treated as a descriptive secondary item.

## 3. O8 audit (frozen calibrator; DT5-ML-DESIGN §2)
- Methods on T = √n·D_n, same samples: uniform q^U_n, LF 0.73139, λ-interpolations λ ∈ {0.25, 0.5, 0.75}, learned q_θ. α = 0.05 only.
- Audit nulls (12): W(0.02), W(0.1), pin κ ∈ {0.85, 1.4, 2.5}, ramp, clipjit 0.002 / 0.05, relu0+jitter, normal, beta(2,2), uniform. Alternatives: training alternatives + mix(2.2). n ∈ {200, 800, 3200}; **R = 50,000** (raised from 20,000 after review: at 20,000 an exact-5% cell passes only ~56% of the time); root 2026110671.
- Honest scope: "a frozen, development-informed calibrator evaluated on independently generated audit samples at prespecified (law, n) cells, with simultaneous Monte Carlo size bounds". The audit laws overlap in mechanism with training (W laws, step adversaries) and DEV (normal, beta); this is parameter/interpolation stress testing, NOT mechanism-held-out validation. The audit panel was chosen with DT5 and DEV knowledge.
- **O8-R1 (size, three labels):** one-sided CP-95% bounds, Bonferroni over method × null cell (6 × 36 = 216, so the label of whichever method is reported is covered). Certified within tolerance = every upper bound ≤ 5.5%; demonstrated above tolerance = some lower bound > 5.5%; otherwise inconclusive. "Not certified" is not "oversized".
- **O8-R2 (deliverable):** (a) learned certified → deliverable (2), with its paired power gain over the most powerful CERTIFIED constant baseline (LF = λ 1, or a certified λ < 1), not only over LF (a gain over LF is guaranteed by q_θ ≤ q^LF). "Learning helped" only if the mean paired gain over that baseline across the alternatives at n = 800 is ≥ 2 pt. (b) learned not certified, some λ certified → report the smallest certified λ as a stage-1 baseline (no further tuning on audit data). (c) only LF certified → deliverable (4), negative for learned calibration at this budget. (d) LF not certified or nothing certified → report as is; no deliverable. (e) uniform certified → report it as a contradiction of DT5 to be explained. No new model or λ tuning on audit results without a new independent audit.

## 4. Predictions (with what informed them)
Adjudication: a prediction "p ≥ a" is SUPPORTED if the CP-95% lower bound > a, CONTRADICTED if the upper bound < a, else INCONCLUSIVE (reverse for "p ≤ a"). A compound prediction is supported only if every component in every listed cell is supported. Paired differences use the Newcombe interval.
- **PM1** (informed by DEV quick §6): UniDip α = 0.05, on N1, N2(σ ≤ 0.01), all N2-pin, at both n (16 cells), three components scored separately: (a) uniform arm P(k̂ > 1) ≥ 10%; (b) LF arm ≤ 5%; (c) M-R1 declared.
- **PM2** (theory: Thm 4, smooth/log-concave): R0(σ), G: both arms ≤ 2%; M-R1 not declared.
- **PM3** (DEV): U, UniDip α = 0.05: uniform arm > 5%, LF arm ≤ 1%. (DEV suggests pipeline over-splitting above the per-test 5%; the root and call logs will show the mechanism. Not assumed.)
- **PM4** (o5_hx diagnostic): DipMeans on nulls: both arms ≤ 1%. Weak by design (nominal 0.1%; DEV found no power on alternatives either); reported with its price-law rates.
- **PM5** (DT5 price-of-validity): UniDip α = 0.05 on mix(3), mix(4), tanh at n = 2000: LF arm ≥ 80%; on mix(2.5) at n = 500 the paired difference uniform − LF ≥ 10 pt. Laws as coded in `o5_ml.sample_1d`: mix(s) = N(0,1) + s·Bernoulli(½); tanh = tanh(3·U(−1,1)); clusters = non-noise labels.
- **PG1** (gap selector, not a dip finding; inconclusive allowed at R = 200): gap P(k̂ > 1) ≥ 20% on N1, N2(0.002), A_N2 at n = 500 (DEV partial: 45% / — / —); ≤ 10% on U, G, logistic.
- **PO1**: LF passes O8-R1 (DT5: LF ≤ 4.5% on similar laws at n ≤ 12,800).
- **PO2**: uniform and λ = 0.25 are "demonstrated above tolerance" (DT5: uniform 26-30%).
- **PO3** (uncertain; training log at the last iteration shows adversarial nulls 5.5-6.6% at n = 800 on batches of 256, a noisy basis): learned is NOT certified within tolerance. A lean, not a strong prediction, and not a forecast that its true size exceeds 5.5%.
- **PO4**: if learned is certified, its mean paired gain over the most powerful certified constant baseline at n = 800 is < 2 pt ("learning did not help").

## 5. Interpretation limits
- O5-ML measures what changes when only the p-value calibration is swapped inside fixed pipelines; it does not say the pipelines are valid or invalid as a whole.
- The audit set is finite; passing O8-R1 is not a validity theorem.

## 6. DEV runs seen before this file (disclosure)
Added (~15:43Z): o5_ml --quick rerun 15:43Z on the final 19-law list (roots 0662/63) to test `code/o5_analysis.py` (K = 70, 19 declared in DEV); o8_audit --quick 15:42Z (root 0673, quick calibrator) to test the three-label rule.

o5_smoke 11:25Z (root 0662); o5_ml --quick 11:33-11:35Z (roots 0662/63, 60 reps; mix2 bug fixed after); o5_hx 11:41Z (root 0664); o5_gap --quick 11:50-12:11Z partial (root 0666, 3 cells: N1 n100 32.5%, n500 45%, N2(0.002) n100 52.5%); o8_calibrator --quick 11:38Z (0672) and full training 11:40-12:15Z (0670); o8_audit --quick 11:39Z (0673, on the QUICK calibrator, not the frozen one). DT5 CONFIRM results seen 15:35Z.
