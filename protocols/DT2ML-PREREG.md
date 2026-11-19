# DT2ML-PREREG — Step 4b: the learned adversary (O7). DEV tier for the search; fresh-seed re-estimates for the read-out.

Frozen 2026-11-05 (file last modified 19:47:21 UTC), before any 4b output (result file last modified 19:48:41 UTC).
Roots (declared in CV-6):
- **2026110620** = training draws (torch generator seed + numpy for the dips);
- **2026110621** = fresh re-estimation of the frozen learned laws;
- **2026110622** = the non-ML step-family baseline.

## 1. Target
The finite-n worst case of the uniform-calibrated 5% dip test at **n = 200** (and a check at n = 1000), using the own-n uniform quantile q̂ from `deposits/dt1b/ref_n{n}.npy` (R = 1e6).

## 2. Generator (shape-constrained; every law is unimodal by construction)
- A density on [0, 1] that is piecewise constant on K = 64 equal bins. The mode bin m is fixed per run.
- Log-heights: ℓ_j = ℓ_m − Σ_{i=j}^{m−1} s_i for j < m (rising to the mode), and ℓ_j = ℓ_m − Σ_{i=m+1}^{j} t_i for j > m (falling), with s_i, t_i = softplus(θ) + δ ≥ δ.
- **Unrestricted track:** δ = 0 (flat stretches allowed; positive control: it should rediscover a flat shoulder + pile-up).
- **Curved track:** δ = 0.05 (each step changes the density by ≥ 5%; no flat stretch longer than one bin). This tests whether the excess needs exact or near flatness.
- Mode positions: m ∈ {63 (right edge), 48 (x ≈ 0.75), 32 (centre)}.
- Inits: (i) uniform (θ ≈ −6); (ii) a peaked start (θ = 0).

## 3. Training signal
- **Score-function (REINFORCE) gradient** of J(θ) = E_θ[σ((√n·D − √n·q̂)/τ)], τ = 0.05 (a smoothed rejection indicator).
  - The log-likelihood of a dataset is Σ_k c_k·log h_k − n·log Z (c_k = bin counts), so ∇J = E[(s − b)·∇ log-lik], with b the batch mean.
- D is computed by diptest (CV-1).
- Adam, lr 0.05; 300 iterations; batch 2,048 datasets per iteration.
- The reported training objective is the batch rejection rate (indicator) at each iteration.

## 4. Non-ML baseline (same budget order)
A grid over the two-block family (shoulder + edge pile-up, mode at the right), with w ∈ 30 log-spaced values in [0.01, 0.5] and ε ∈ 30 log-spaced values in [0.001, 0.2] (ε < w), at n = 200 and 4,000 datasets per cell.

## 5. Read-out (fresh seeds; target: a reproducible extremal candidate plus an explanation)
- Freeze the best law of each run (the final θ).
- Re-estimate its rejection at n = 200 and n = 1000 with 100,000 datasets on root 2026110621.
- Also re-estimate the best baseline law the same way.
- Save each learned density (64 heights) and a description: the block masses, the largest near-flat stretch, and the pile-up size.

## 6. Questions this answers (feeds the kill table)
- **Q-a:** does the unrestricted adversary rediscover the flat-shoulder + pile-up (K1/K4)?
- **Q-b:** does it beat the best 2-block law at n = 200, i.e. is the finite-n worst case outside the 2-block family?
- **Q-c:** with the mode at the centre, does it find any law > 5% (K9)?
- **Q-d:** does the curved track find anything > 5%?

## 7. Recorded S predictions
- Q-a: yes (a long flat shoulder with an edge pile-up of mass ≈ 0.1-0.2).
- Q-b: the learned best is within ±1 point of the best 2-block law, about 26-28% at n = 200.
- Q-c: with the mode at the centre, it moves the mass to one side (an unequal interior plateau): about 15-20%.
- Q-d: curved track ≤ 5%, or only slightly above (at most 6-8% at n = 200, from near-flat bins).
