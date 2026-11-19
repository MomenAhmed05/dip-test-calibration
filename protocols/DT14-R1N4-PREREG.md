# DT14-R1N4-PREREG v2 — Step 15 WP-D at n = 4: a certified valid dip cutoff (R1)

Status: v2 FROZEN 2026-11-08 17:31Z (`date -u` 17:30:59 at hashing). The v1 draft was revised after one review round; v2 contains all its fixes.

## 1. What will be certified
For the three prespecified cutoff pairs (α, d_α) = (10%, 0.2335), (5%, 0.2420), (1%, 0.2488) we will conduct a simultaneous Monte Carlo
certification of the assertions
  sup_{F ∈ U_c} P_F(D_4 >= d_α) <= α,
where U_c = continuous unimodal laws in the Hartigan sense (CDF convex before a mode, concave after it). The certification procedure has
false-certification probability at most 0.001 under the stated iid simulation model (repeated-simulation coverage; NOT a posterior
probability and NOT an unconditional deterministic proof). Inputs (§3) are conditional.

## 2. Why it matters (R1 at n = 4)
- Uniform-calibrated cutoffs (dip >= 0.1875 / 0.2072 / 0.2317) are INVALID at n = 4: Thm N4-LB (`DT14-N4.md` §1, PROVEN-written) gives worst size >= 13.746 / 8.596 / 3.264% > α.
- The specified asymptotic-cutoff rule ("LF": reject iff D_n >= c_α/√n, c_α = 0.6684 / 0.7314 / 0.8599) has cutoffs 0.334 / 0.366 / 0.430 at n = 4, above the largest possible dip 1/4: it never rejects at n = 4.
- Power comparison (deterministic, rejection-region inclusion): at n = 4 each proposed rule weakly dominates the specified asymptotic-cutoff rule in power over all alternatives, and strictly dominates it at every alternative G with P_G(D_4 >= d_α) > 0. This is not a claim against every finite-sample calibration.
- If Conjecture N4* holds, the optimal valid cutoffs would be dip 0.2015 / 0.2231 / 0.2440; ours are larger (conservative) because the only proved upper bound is Lemma N4.

## 3. Chain of reasoning (inputs named, tiers)
1. n = 4 dip (H&H 1985: the dip is an increasing function of 1 ∨ T; explicit form derived by us and checked = B0): D_4 = max{1/8, T/(4(T+1))}, T = s2/max(s1, s3). Hence for d ∈ (1/8, 1/4): {D_4 >= d} = {T >= T_d}, T_d = 4d/(1-4d) = 467/33, 121/4, 622/3 (all > 1). B0 check (DEV root 2026110815): 20,000 samples (uniform, normal, multi-scale, near-tie), max |diptest − formula| = 5.6e-17; LP dip (3,000) max diff 2.1e-9 (LP tolerance). Written derivation of the explicit form: B0 item, pending.
2. Lemma N4 v2 (DT13-PREREG §4; PROVEN-written): pathwise {T_F >= r} ⊆ E_u under the quantile coupling, u = F(m) for a fixed population mode m; E_u = {φ(u) >= r}, φ(v) = A/g1 + (g2 − A)/g3, A = min(max(v − b, 0), g2), (a,b,c,d) order statistics of 4 iid U(0,1). Hence sup_F P_F(T >= r) <= sup_{u ∈ [0,1]} P(E_u) = sup_{u ∈ [1/2,1]} P(E_u) (reflection (a,b,c,d) → (1−d,1−c,1−b,1−a), u → 1−u).
3. φ is monotone in v for each sample (constant, affine with slope 1/g1 − 1/g3, constant) ⇒ ∪_{v ∈ [L,R]} E_v = E_L ∪ E_R. Grid u_j = 1/2 + j/200, j = 0..100; blocks B_j = [u_{j−1}, u_j] cover [1/2, 1] including both ends. The simulated indicator is the single union indicator 1{E_L ∪ E_R}.
4. Ties: for continuous F and r > 1, {T = r} ⊆ {s2 = r s1} ∪ {s2 = r s3}, nontrivial linear equalities on each ordering; conditioning on three observations, an atomless fourth hits a single value with probability 0. So P_F(D_4 = d) = 0 for d ∈ (1/8, 1/4); ">=" and ">" give the same rejection probability under every continuous null (not claimed for discrete alternatives).
5. MC: K_kj ~ Binomial(N, p_kj), p_kj = P(E_{u_{j−1}} ∪ E_{u_j}) at cutoff k; one-sided CP upper limit U_kj = Beta^{-1}(1 − η; K + 1, N − K) (= 1 if K = N), η = δ/300, δ = 0.001; union bound ⇒ P_MC(∀k,j: p_kj <= U_kj) >= 0.999 (same quadruples for all events is fine).
6. Numerical safeguard: event evaluated division-free, A g3 + (g2 − A) g1 >= r g1 g3 (1 − 10^{-12}), i.e. counted INCLUSIVELY (more hits can only raise the bound). CP covers sampling error only.

## 4. Run
- Script `code/dt14_r1_n4.py` (sha recorded at freeze). Command: `python code/dt14_r1_n4.py 2026110830 20000000` (from the repository root).
- Root 2026110830 (first Step 15 CONFIRM root; single use; banned afterwards). N = 2·10^7, chunks 2·10^6, numpy PCG64 (numpy 2.3.5 recorded in the output).
- Output records all 300 counts and limits.
- DEV before freeze (roots 2026110811, 0813, 0814, 0816; independent of the CONFIRM root): max block P ≈ 0.0960 / 0.0459 / 0.0067 at [0.500, 0.505]; quick N = 4·10^5 max CP 0.0983 / 0.0477 / 0.0073.

## 5. Decision rule (frozen)
Simultaneous PASS = ∩_k {max_j U_kj <= α_k}; each level also reported separately. FAIL at a level ⇒ no claim at that level; no re-run (a new prereg version would be required: stop rule).

## 6. Predictions (before the run)
max CP ≈ 0.0962 / 0.0461 / 0.0068 (all PASS); argmax block [0.500, 0.505] at all levels.

## 7. What a PASS does and does not mean
Does: valid n = 4 dip tests at each level, conditional on items 1-2 of §3, with false-certification probability <= 0.001; strict power gain over the specified asymptotic-cutoff rule as scoped in §2.
Does not: optimality (depends on N4*, OPEN); anything at n >= 5; a deterministic proof.
