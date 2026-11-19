# DT2b-PREREG — Step 4 farm part 2: the pinning curve and the shrinking jump. DEV tier.

Frozen 2026-11-05 (file last modified 19:40:44 UTC), before any part-2 output (result file last modified 19:45:10 UTC). Root **2026110611** (DEV; declared in CV-6): reference (0, n, chunk), laws (1, f, n, chunk).

## 1. Laws
- **Pinning family P(κ, λ):** pile-up at the right end, width ε = λ/√n and excess mass e = κ/√n, so w = ε + e(1 − ε). Blocks: (0, 1 − ε) with mass 1 − w; (1 − ε, 1) with mass w (modal).
  - κ ∈ {0, .2, .4, .6, .8, 1, 1.2, 1.4, 1.7, 2, 2.5, 3, 4}; λ ∈ {0.1, 1, 4}; n ∈ {3200, 12800}.
  - κ = 0 is the exact uniform law (the anchor).
- **Shrinking jump S(c):** blocks (0, .5) with height h1; (.5, .995) with height h1(1 + c/√n); (.995, 1) with mass .1 (modal). The total mass on (0, .995) is .9.
  - c ∈ {0, 1, 2, 4, 8, 16}; n ∈ {3200, 12800}.
  - Bounds: "merged" = BF with one block of mass .9; "separated" = BF with the two blocks' own masses.

20,000 datasets per cell. Uniform reference R = 1e6 per n.

## 2. Readings fixed in advance (α = 5%)
- **(a) λ-collapse:** at n = 12800, for every κ, max − min over λ ≤ 2.0 points → "COLLAPSE (λ irrelevant at this precision)". Otherwise "λ MATTERS" (list the κ).
- **(b) n-stability (λ = 1):** |g(3200) − g(12800)| ≤ 2.0 points for every κ → "LOCAL LIMIT STABLE".
- **(c) No overshoot:** max over all P cells ≤ r_sup(n) + 1.0 point, where r_sup(n) = 1 − F_M(2·√n·q̂_n) is the global BF sup at the own-n quantile → "NO OVERSHOOT".
- **(d) Shrinking jump:**
  - every cell lies within [separated − 1, merged + 1] points;
  - the rejection is nonincreasing in c (allowing 1 point of noise);
  - |3200 − 12800| ≤ 2 points → "SANDWICH + LOCAL LIMIT".

## 3. Recorded S predictions
- (a) COLLAPSE.
- g(κ) at n = 12800: κ = 0 → 5%; κ = 0.85 → ≈ 19%; κ = 1.7 → ≈ 30%; κ = 4 → ≈ 31.5%. Monotone.
- (b) STABLE.
- (c) NO OVERSHOOT.
- (d) c = 0 → ≈ 24.8% (merged); c = 1 → ≈ 15%; c = 2 → ≈ 10%; c = 4 → ≈ 6%; c = 8 → ≈ 4.5%; c = 16 → ≈ 4%. SANDWICH + LOCAL LIMIT.
