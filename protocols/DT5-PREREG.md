# DT5-PREREG — Step 7 size and power study of dip calibrations

**Status: FROZEN** ~03:34Z on 2026-11-06 (`date -u` 03:33:01 when the code hashes below were taken). Draft history: started 2026-11-06 00:59Z; revised after one review round.
Design type: **pilot-informed, frozen before the independent main simulation**. Two DEV quick runs (roots 2026110652/53, 400 samples per cell, 03:21Z and 03:31Z) were seen; §8 lists what they showed and which parts of this file were written after them.

Pipeline (frozen): `code/dt5_power.py` (simulator, sha 1d0f2c52), `code/dt5_analysis.py` (intervals + rules, sha b4fbe8ad), `code/hartigan2000.py` (sha a0d21e73), `code/dt0c_screen.py` (`dips_of`, sha 7fe20a4f). Versions: Python 3.12.10, numpy 2.3.5, scipy 1.17.1, diptest 0.12.0.

## 1. Methods (same datasets for all; paired)
Primary (one statistic, the dip D_n, at four critical values):
1. **uniform** — own-n uniform reference quantile. The procedure under repair; reported as a reference, never ranked. Its advantage over LF is called a "uniform-calibrated rejection advantage", not a power gain.
2. **LF** — c_α/√n, c_α = ½F_M⁻¹(1−α) = 0.73139526 / 0.85987427 (Balabdaoui-Filali Table 5 halves; DT4-CERT). Scope: sharp, pointwise asymptotic level over finite-step unimodal laws (Theorem 3 + Prop. 3). Not proved: finite-n validity, validity over all unimodal laws, uniformity.
3. **kuiper2** — ½ × a conservative order-statistic upper bound (one-sided 1 − γ, γ = 10⁻³ per (n, α); total calibration failure budget ≤ 0.018 over 18 cutoffs) on the own-n Kuiper quantile from 10⁶ uniform samples. Scope: finite-n valid for continuous unimodal F (Davies & Kovac 2004 §3.6: order-1 Kuiper distance to unimodal = 2 × dip). The point quantile is also recorded (`kuiper2_point`).
4. **dkw** — √(log(2/α)/(2n)) (DKW-Massart). Finite-n valid, very conservative.

Secondary diagnostic, appendix only: **sse_asymptotic** and **sse_uniformcal** = "an experimental dip-selected, single-shoulder excursion procedure inspired by Hartigan (2000), not an implementation of his full antimode procedure" (T = N·D_s/√m, m ≥ 5 else no rejection; Kennedy = Kuiper reference 1.74725995 / 2.00091812, or own-n uniform critical values from 200,000 samples). Neither reference is asserted to give composite-null size control. Any failure concerns this implementation.

## 2. Laws (37; frozen manifest written into `power.json`)
- **Nulls (9):** uniform; W(0.1, .005); N(0,1); Beta(2,2); 0.9U + 0.1Beta(0.5,1); W(0.02, .001); pin(1.4, .1) and pin(2.5, 1) (pinning sequences: top width λ/√n, excess k/√n; the law changes with n); smoothW(0.02, .001, 10⁻⁴) = W(0.02, .001) + N(0, 10⁻⁸) noise (unimodal by a one-sign-change argument for f′; finite-n interpretation only — not claimed to be asymptotically least favourable).
- **Unimodal controls** (counted as nulls in R1): ½N(∓δ,1) with δ ∈ {0.8, 1}; asymmetric mixtures with 1 mode; valley step laws with ρ ≥ min(1, s).
- **Alternatives (bimodal):** symmetric δ ∈ {1.05, 1.1, 1.25, 1.5, 2, 2.5, 3}; asymmetric w N(−δ,1) + (1−w) N(δ,1), w ∈ {0.1, 0.25}, δ ∈ {1.25, …, 3} with 2 modes; valley laws with ρ < min(1, s); satellite 0.9 N(0,1) + 0.1 N(3, 0.25²).
- Mode counts: closed form for symmetric and valley laws; for Gaussian mixtures by sign changes of f′ on a 2·10⁶-point grid with a near-tangency flag (numerical, labelled; no law was flagged).

## 3. Design
- n ∈ {50, 100, 200, 400, 800, 1600, 3200, 6400, 12800}; α ∈ {0.05, 0.01}.
- Reference: 10⁶ uniform samples per n (dip + Kuiper; shoulder T on the first 200,000). Root **2026110650**, spawn (0, n, chunk).
- Laws: **50,000** samples per (law, n) (raised from 20,000 after review: at 20,000 a method of exact size 5% fails the R1 certificate with probability 12%). Root **2026110651**, spawn (1, law_index, n, chunk).
- Expected run time ≈ 1-1.5 h on 11 workers. Per-job parts are saved; a restart skips finished jobs.
- Resource fallback (fixed now, independent of results): if the run is interrupted and cannot be resumed within 24 h, the shoulder statistic on laws at n ≥ 6,400 is dropped first; replication counts and laws are never changed because of results.
- Numerical inputs (CV-11): predictions use F_M at x ≥ 1.09, inside the published table.

## 4. Readouts (all by `dt5_analysis.py`)
Per (law, n, α, method): count, rate, two-sided Clopper-Pearson 95% interval (conditional on the realised reference cutoffs). Paired: LF − kuiper2 extra rejections and uniform − LF extra rejections (nested ⇒ one discordance probability, CP interval). Power: first registered n with point power ≥ 50% / 80%, and first n whose CP lower bound ≥ 50% / 80% ("not reached" allowed; no monotonicity assumed). Shoulder diagnostic: failed rows, eligible fraction, median m, median R = D_s/dip.

## 5. Decision rules (frozen)
- **R1 (size tolerance on the panel).** Bounds are the endpoints of two-sided CP-95% intervals. For method M at (n, α), over all null and control laws: **certified within tolerance** if every upper bound ≤ α + Δ_α (Δ = 0.5 pt at 5%, 0.25 pt at 1%); **demonstrated above tolerance** if some lower bound > α + Δ_α; else **inconclusive**. Failure to certify is not evidence of inflation. The laws whose lower bound exceeds α are listed separately. Shoulder rows with any failed sample are INVALID for that cell.
- **R2 (persistent finite-sample overshoot of LF on pinning).** Four claims (2 sequences × 2 levels). Per claim: one-sided binomial p-values of "rate > α" at n = 6,400 and 12,800; claim p = the larger; Holm across the 4 claims; alarm if Holm p < 0.05. An alarm is finite-sample evidence only; it does not disprove asymptotic validity (C4), and c_α is never inflated in response (DT5-DESIGN §2).
- **R3 (comparisons).** Equal-level power comparisons only between methods certified at that (n, α). Methods are labelled by scope: theorem-backed (LF: asymptotic, step laws; kuiper2, dkw: finite-n), or no validity theorem (shoulder diagnostics). All curves are shown.
- **R4 (sanity).** The realised order of cutoffs LF < kuiper2 < dkw is checked per (n, α); whenever it holds, rejections must be nested dataset by dataset. A nesting violation = implementation error, the run is void until explained. An order violation triggers a numerical investigation.
- **R5 (headline, replaces "most powerful valid").** Report LF's paired gains over kuiper2 and dkw with their intervals, LF's R1 status, and its theoretical scope; and the uniform-calibrated rejection advantage over LF as the price of validity. No optimality claim among unimodality tests.

## 6. Predictions (S, before the main run)
- P1: LF certified within tolerance at all n ≥ 800 (theory: W(0.1) limit 3.1% / 0.5%; W(0.02) 4.6% / 0.9%; smooth nulls ≈ 0). The DEV quick run (§8) leaves this uncertain at n = 12,800 for the W(0.02)-type nulls.
- P2: uniform demonstrated above tolerance at every n ≥ 100 (W-type and pinning nulls ≈ 24-31% at 5%).
- P3: R2 no alarm.
- P4: kuiper2 and dkw certified everywhere.
- P5: both shoulder diagnostics demonstrated above tolerance at n ≥ 400 (informed by DEV, §8).
- P6: kuiper2 needs about 1.4× as many observations as LF for 80% power on symmetric mixtures (heuristic (0.871/0.731)²; LF needs about 29% fewer).
- P7: R4 passes (0 nesting violations).

## 7. Interpretation limits (fixed now)
- Size results hold for this panel only; they do not prove validity elsewhere.
- The shoulder rows test OUR one-shoulder simplification. Wording if they fail: "a dip-selected single-shoulder simplification of Hartigan's (2000) statistic did not control size on smooth unimodal laws in our simulations".
- The uniform-calibrated rejection rates on alternatives are not power at level α.

## 8. What the DEV quick runs showed and what was written after them
- Run 1 (03:21Z, older code): shoulder statistic at 5% on N(0,1) 5.2 / 7.0 / 9.8 / 15.2% at n = 50 / 400 / 3,200 / 12,800; R = D_s/dip = 2.00 in every null cell; LF 2.5-5.0% on W-type and pinning nulls at n = 12,800.
- Run 2 (03:31Z, final code): R1 status — uniform ABOVE at all n; LF certified at n = 50, 200, 800, inconclusive elsewhere (max null rate 6.0% at n = 12,800, smoothW, 400 samples, SE ≈ 1.2 pt); kuiper2/dkw certified or inconclusive; shoulder diagnostics ABOVE from n = 100. R2 no alarm. R4 18/18, 0 violations.
- Written before any DEV run: methods 1-4, the R1 5% margin and the R4 order idea (00:59Z draft). Written after run 1: P5, the R2 threshold. Written after run 2 and the review round: the R1 three-label form and 1% margin, R2 Holm form, R3 scope labels, R5, the 50,000 replication count, the two added laws, the conservative Kuiper cutoff, the failed-row policy.
