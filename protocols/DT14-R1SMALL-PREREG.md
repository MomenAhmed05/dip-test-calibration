# DT14-R1SMALL-PREREG v1 — Step 15 WP-D (T15-D): certified dip cutoffs below the fix at every n = 8..49 (R1)

Status: v1 FROZEN 2026-11-08 19:56Z (`date -u` 19:56:13 at hashing), after one review round whose fixes are listed in §8.

## 1. What will be certified
For every n ∈ {8, 9, …, 49} and α ∈ {10%, 5%, 1%}, a dip cutoff d_{n,α} chosen by a frozen rule (§3) such that

  sup_{F ∈ U_c} P_F(D_n ≥ d_{n,α}) ≤ α,

U_c = continuous unimodal laws (Hartigan sense: CDF convex before a mode, concave after it), D_n = Hartigan dip.
Guarantee: P_CONFIRM(some reported PASS cell has sup_F P_F(D_n ≥ d) > α) ≤ δ = 0.001, under the iid uniform simulation
model and conditional on §2 items 1-2. This is repeated-simulation coverage, not a posterior probability and not a
deterministic proof. Only the 42 × 3 listed cells are covered; nothing is interpolated to other n.

## 2. Chain (inputs and tiers)
1. Cor U1′ (DT7-UNIF §1, PROVEN): samplewise D(F_n) ≤ S̄_{n,I_j} for F ∈ U_c with u = F(mode) ∈ I_j
   after reflection (a distribution-level reflection, Step B, not sample-chosen), so sup_F P_F(D_n ≥ d) ≤ max_j P(S̄_{n,I_j} ≥ d).
   15 closed blocks with edges 0.5, .55, .6, .65, .7, .75, .8, .85, .9, .93, .95, .97, .98, .99, .995, 1 (cover [½, 1]; shared
   endpoints belong to both blocks). Pieces [a, b] = I_j ∩ [x_(k), x_(k+1)] with x_(0) = 0, x_(n+1) = 1, closed, singletons kept.
2. Numerical certificate (module docstring of `code/dt9_block_cert.py`): piece_ub = band error of the returned
   knot vector + PROVED shape-repair cost (h = Σδ_k(t_k − s)⁺ − Σδ′_k(s − t_k)⁺ makes g + h ∈ C_[a,b], |h| ≤ Σδ + Σδ′) +
   floating-point rounding allowance 16u(Σ|slopes| + K·rep + 1) + 1e-9. This is a floating-point certificate with an explicit
   rounding allowance, NOT interval arithmetic (stated as such; same tier as DT9). LP failure ⇒ +∞ (counts as exceedance).
3. Branch and bound per block (`code/dt14_r1_small.py`): a node is pruned only when its certified bound is ≤ max(FLOOR, best
   leaf so far); node bounds dominate every piece inside (C_[a_i,b_j] ⊆ C_piece). So the stored R_j ≥ S̄_{n,I_j} whenever
   R_j > FLOOR, and FLOOR = 0.999 × (lowest ladder cutoff) < every ladder cutoff − 1e-12. DEV check (root 2026110821, n = 8/20/49,
   150 samples each): R_j ≥ (full piece enumeration) − 1e-12 in all cells; max difference 3e-13.
4. Event: exceedance counted as R_j ≥ d − 1e-12 (inclusive; can only raise counts). The rejection rule is "reject iff D_n ≥ d".

## 3. Frozen cutoff rule (fixed sequence)
- Ladder per (n, α): d_l = ρ_l · c_α^lo / √n, ρ_l = 0.900, 0.905, …, 0.995 (20 values), c^lo = 0.667 / 0.7313 / 0.8585 = LOWER
  ends of the certified LF enclosures (so every ladder cutoff is strictly below the true LF constant's cutoff).
- Counts K_lj on the same N CONFIRM samples; U_lj = one-sided Clopper-Pearson upper limit at level η (BetaInv(1 − η; K + 1, N − K),
  1 if K = N), η = δ / (42 · 3) = 7.94e-6. PASS at l iff max_j U_lj ≤ α. Reported cutoff = the FIRST passing l; none ⇒ no claim.
- Why no factor for blocks or ladder (checked by S): counts, hence U_lj, are nonincreasing in d samplewise; the true
  p(d) = max_j P(S̄_j ≥ d) is nonincreasing; let l* be the largest ladder index with p(d_{l*}) > α (deterministic). A wrong report
  needs a pass at some l ≤ l*, hence a pass at l*, hence U_{l*,j*} ≤ α < p_{j*}(d_{l*}) at the FIXED block j* attaining p > α:
  a CP non-coverage event of probability ≤ η. Union over 126 cells ⇒ ≤ δ.

## 4. Run
- Root 2026110831 (single use; banned afterwards). N = 100 000 per n (20 chunks × 5 000; chunk k of n uses
  SeedSequence([2026110831, n, k]), numpy PCG64). Script `code/dt14_r1_small.py` (sha recorded at freeze), 11 parallel workers over n.
- Analysis `code/dt14_r1_small_analysis.py 2026110831 100000 deposits/dt14/r1small_2026110831.json` (sha recorded at freeze).
- No re-run, no retuning, no added samples after the first look. A crash is resumed with the same root (existing chunks skipped).

## 5. Decision and reporting
Per cell: certified cutoff or "none". Report every cell (including none and cells whose cutoff is ≥ the explicit two-cluster dip,
"cannot reject": 1/4 for even n, (n − 1)/(4n) for odd n). A "none" cell is a FAIL for that cell (fallback: no claim; the
never-reject test is the only safe fallback). Global result = list of certified cells.

## 6. Predictions (DEV pilot root 2026110819, N = 20 000 at n = 8, 10, 12, 15, 20, 30, 49; 1% steps)
Certified ρ ≈ 0.94-0.95 (n = 8-12), 0.96-0.97 (n = 15-20), 0.97-0.98 (n = 30-49) at all three α; the binding block is the edge
[.995, 1] (block probabilities increase smoothly in u, consistent with C-edge). Cannot-reject cells expected at 1% for n = 8, 9, 10
(and possibly 11). With the finer 0.005 ladder, expect ρ about 0.005 lower than the 1%-step DEV values. All 126 cells expected PASS
at some ρ ≤ 0.99.

## 7. What a PASS means and does not mean
- Means: at each certified (n, α) the dip test "reject iff D_n ≥ d_{n,α}" has finite-sample size ≤ α over all of U_c (with the §1
  guarantee). Because d_{n,α} < c_α^lo/√n ≤ the LF cutoff, its rejection region contains the LF rule's, so its power is ≥ the LF
  rule's at EVERY alternative, strictly where P(d ≤ D_n < LF cutoff) > 0 (deterministic inclusion; no alternative simulation).
  By-product: in the same cells the LF rule itself is certified finite-sample valid (it rejects on a subset).
- Does NOT mean: the cutoff is the smallest valid one (U1′ is loose at small n; the uniform quantile is far lower and invalid);
  dominance over a finite-n optimally calibrated dip test; any power comparison with Dümbgen-Walther 2008 / Rufibach-Walther 2010
  multiscale tests (no rejection-region inclusion; not run here). Wording: "improved certified calibration of the dip statistic",
  not "a more powerful modality test".

## 8. Changes made after review (labels refer to the review questions)
Q2: fixed-sequence ladder frozen, pass = all-block CP ≤ α. Q4: edge block kept as a block (no endpoint replacement); refinement of
the edge was checked in DEV (blocks [.99, .995] and [.995, 1] are within 0.001 of each other, so splitting cannot help); stratified
sampling not used. Q5.1: certificate tier stated (§2.2). Q5.2: lower LF ends used. Q5.3: all 42 integers n = 8..49; odd-n dip bound
stated as an explicit example (can-reject), not a supremum. Q5.4: N chosen for pass probability (DEV: CP ≈ p̂ + 4.5√(p̂/N)).
Q5.5: guarantee wording as §1.
