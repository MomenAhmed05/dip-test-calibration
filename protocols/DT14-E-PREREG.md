# DT14-E-PREREG v1.1 — Step 15 WP-E (T15-E): certified dip cutoffs vs global rules derived from DW08 / RW10

Status: v1 FROZEN 2026-11-09 01:43Z (`date -u` 01:42:30 at hashing), after one review round. **v1.1 AMENDMENT (2026-11-09 01:48Z, `date -u` 01:47:39 at hashing), before any power sample was drawn:** the v1 calibration code capped the RW10 block-criterion search at
ã ≤ 1 (`hi <- min(1, α(A+1)²)`), which is NOT in the RW10 definition or in §3 below (only ã/(A+l)² ≤ 1 is needed); with few
blocks it held the RW tests at ~2% union exceedance for α = 10% and 5% (unfair to RW10). Fixed (`hi <- α(A+1)²`); calibration
re-run on the fresh root 2026110834; the root-0832 deposit is kept in `deposits/defects/mh_calib_root0832_capbug/` and 0832 is
banned. DW08 cutoffs were unaffected by the bug. No other change. The power root 2026110833 was unused at the time of the fix.
This is a **preregistered final simulation following exploratory development**
(DEV roots 2026110825-29 shaped the code and the panel; see §8).

## 1. Scope
We compare rejection probabilities of class-valid dip procedures and specified global unimodality rules derived from DW08 / RW10
simultaneous directional procedures, at 11 prespecified sample sizes between 8 and 49 (N ∈ {8, 10, 12, 15, 20, 25, 30, 35, 40, 45,
49}) and three nominal levels (α = 0.10, 0.05, 0.01). Claims concern this alternative panel, with calibration guarantees and
certification uncertainty stated separately. The comparison does not assess the procedures' localisation performance or establish
dominance over unimodality tests generally.

## 2. Tests (all evaluated on the SAME samples)
| id | rule | guarantee (stated separately from power) |
|---|---|---|
| R1 | D_N ≥ d_cert(N, α), the first passing ladder cutoff of T15-D (`deposits/dt14/r1small_2026110831.json`) | sup over continuous unimodal laws of P(D_N ≥ d) ≤ α, with SIMULTANEOUS certification confidence ≥ 0.999 over all 126 T15-D cells (conditional on Cor U1′ + float certificate). Cells without a certificate: "no test" (never counted as 0). |
| LF | D_N ≥ c_α^hi/√N, c^hi = 0.6685 / 0.7315 / 0.861 | asymptotic (C2/C13). Finite-N valid exactly in the cells where R1 is certified, because d_cert ≤ c^lo/√N ≤ c^hi/√N makes LF's rejection set a subset of R1's. |
| TAB | diptest 0.77.2 dip.test(x)$p.value ≤ α (uniform-table interpolation) | NOT valid over continuous unimodal laws (our C24: size enclosures above α at every listed (n, α), e.g. n = 20: [16.06, 45.39]% at 5%; C6: W(0.1, 0.005) above nominal at every n 4-1000). Reference only, shown separately. |
| DWa | DW08 all intervals, additive correction Γ_jk = √(2(1 + log((N − 1)/(k − j)))) (as coded in modehunt::modeHunting) | finite-sample (DW08 §2) with a frozen-cutoff tolerance calibration (§3) |
| DWn | DW08 all intervals, no Γ | as DWa |
| RWa | RW10 block criterion (A = 10) on ALL intervals, blocks by point count k − j + 1 in classes [2^c, 2^{c+1} − 1] (c = 1, 2, …), indexed LARGEST class first | RW10 block criterion (RW10 §4 notes it applies to all intervals); co-primary; labelled "prespecified all-interval small-sample partition" |
| RWd | RW10 block procedure with modehunt::modeHuntingBlock's interval set (d0 = 2, m0 = 10, fm = 2), blocks rows 1..n.blocks of modehunt::blocks (largest first) | documents default behaviour; no interval exists for N < 20 and two ordered intervals are impossible for N < 23 (structural zeros, reported separately) |
Global rule for DW/RW: with intervals (X_(j), X_(k)) between order statistics of the raw
sample (unknown support, no added endpoints), reject iff some significant-decrease interval (Z_jk < −c) and some
significant-increase interval (Z_jk > c) satisfy X_(k_dec) ≤ X_(j_inc). Implementation: `code/r/mh_stats.R` (vectorised; decisions
equal modehunt's on 300/300 DW08 and 223/223 RW10 DEV inputs; modehunt::modeHuntingBlock itself raises "incorrect number of
dimensions" on 77/300 small-N inputs — our port never errors; an error would never be read as a non-rejection).

## 3. Calibration (`code/r/mh_calib.R`, root 2026110834 (v1.1; was 0832 in v1), M = 10⁵ uniform N-samples per N, the SAME N-point construction)
- DELTA_C = 0.001/(11 · 3 · 4) per frozen cutoff (Bonferroni over N, α and the 4 calibrated methods).
- DWa / DWn: tolerance order statistic T_(k), k = smallest index with BetaCDF(1 − α; k, M − k + 1) ≤ DELTA_C: the frozen
  cutoff has uniform exceedance ≤ α with calibration confidence ≥ 1 − DELTA_C.
- RWa / RWd: largest ã with empirical union exceedance ≤ α on half 1 (thresholds q_l = type-1 empirical quantiles of the block maxima
  at level 1 − ã/(A + l)²); certification gate on half 2: frozen ladder ã·1.1^{−j}, j = 0..80, first step whose one-sided
  CP upper limit (level DELTA_C) of the half-2 union exceedance is ≤ α (fixed sequence; union probability nonincreasing). No pass ⇒
  the test never rejects (recorded).
- Power intervals (§5) are conditional on the frozen cutoffs; calibration uncertainty is the stated DELTA_C budget.

## 4. Laws (`code/r/mh_power.R`, LAWS; executable definitions; nulls are diagnostics only)
Nulls: U(0,1); N(0,1); Exp(1); t₃; Beta(2,5); trapezoid U(−1,1) + U(−0.3,0.3); W(0.1, 0.005) = with prob 0.1 uniform on
[1 − 0.005, 1], else uniform on [0, 1 − 0.005].
Alternatives: p N(0,1) + (1−p) N(Δ, s²) with (p, Δ, s) = (0.5, 2.2, 1), (0.5, 3, 1), (0.5, 4, 1), (0.5, 5, 1), (0.5, 8, 1), (0.3, 3.5, 1),
(0.3, 5, 1), (0.5, 3, 0.5); ½U(0,1) + ½U(1+g, 2+g), g = 0.5, 1, 3; ½N(0, 0.05²) + ½N(1, 0.05²); local valley on [0,1] with density ∝
1 − h·max(0, 1 − |x − ½|/w), (w, h) = (0.10, 0.6), (0.10, 1.0), (0.25, 0.6), (0.25, 1.0) (multiscale-friendly family
chosen by geometry, not by pilot ranking). Modality documentation in the result: population mode count and locations by grid root
finding (grid 10⁻⁴ of the range), valley depth = 1 − f(antimode)/min(f(modes)), population dip approximated by the dip of the
quantile sample F⁻¹((i − ½)/10⁵) (numerical, labelled).

## 5. Run and analysis
- Power root 2026110833; set.seed(root + 1000 N + law index); N_sim = 10⁴ samples per (law, N); 11 workers (one per N).
- Output per cell: counts per test and α; the 2×2 paired table (R1 vs each test) per α.
- Rates: two-sided 99% Clopper-Pearson (exact; zero counts fine). Paired differences R1 − X: Newcombe (1998) method 10 99% interval
  (approximate). All intervals POINTWISE DESCRIPTIVE: no family-wide superiority claim. Primary level for the
  headline table: α = 0.05; the other levels in full tables.
- Reporting: certificate availability table first; R1-vs-X contrasts only on cells where R1 is certified; structural zeros (RWd,
  N < 23), estimated zeros and "no test" cells labelled differently; TAB shown in a separate panel.
- Mathematical facts vs empirical hypotheses: R1 ≥ LF in every certified cell is BUILT IN (nested rejection sets);
  everything else in §6 is a hypothesis.

## 6. Hypotheses (from the DEV pilot, N_sim = 200-300; NOT promises)
- N ≤ 20: all valid tests ≈ 0-2% except tight clusters; TAB much higher (and invalid).
- N = 49: R1 − LF ∈ [0, 4] points; R1 ≫ DWa, DWn, RWa, RWd on separated laws; near-transition and valley laws low for all valid tests.

## 7. Not covered
Critical-bandwidth / excess-mass bootstrap tests; localisation quality of DW/RW; N outside the 11 values.

## 8. Development disclosure
DEV runs (roots 2026110825-29) informed: the global rule check, the ports, the decision to add tight clusters, N(0,1)+N(8,1), U-gap 3
(after the N = 20 pilot showed near-zero power everywhere) and the valley family. Fresh roots 0834 (calibration, v1.1) and
0833 (power) are used once; 0832 was consumed by the defective v1 calibration (defect record kept).
