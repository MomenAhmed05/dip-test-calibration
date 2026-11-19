# DT0c-PREREG — Step 2b early breadth screen (K10 + realistic pile-up shapes). DEV tier only.

Frozen 2026-11-05 (file last modified 17:02:48 UTC), before any screen output existed (result file last modified 18:27:38 UTC). Conventions: CONVENTIONS.md v1 (CV-1..CV-10).

**Tier:** DEV. These numbers can never confirm a claim. The laws chosen for confirmation are re-run later on fresh roots, with the code paths re-audited.

## 1. Statistic, calibration, estimand
- D_alg (CV-1), calibrated by our own uniform reference H_n with R = 1,000,000 per n (CV-3), root 2026110540, spawn (0, n).
- α ∈ {0.05, 0.01}. Reject iff D ≥ c (CV-4).
- Bracket: γ = 0.001 over the 8 (n, α) pairs, i.e. 1.25e-4 per pair, split equally between the two tails.
- Each cell: 20,000 datasets, root 2026110540, spawn (1, f, n). Every dip is stored (`deposits/dt0c/dips_f{f}_n{n}.npy`).
- Reported per cell: the point rejection at q̂; CP-lower at X_(u); CP-upper at X_(l).

## 2. Grid
n ∈ {100, 200, 500, 1000}. Families f = 0..48 (§3). 49 × 4 = 196 cells.

## 3. Families, roles and membership certificates
Notation: Z, E ~ N(0,1) independent; IQR_Z = 1.3489795; B ~ Beta(a, 1) with F_B(x) = x^a; q_a = 0.9^{1/a}; IQR_{B,a} = 0.75^{1/a} − 0.25^{1/a}; IQR_Exp = ln 3. "+ k% noise" = + σE with σ = k/100 × the IQR of the uncensored component.

| f | family | role | membership (CV-10 rule) |
|---|---|---|---|
| 0 | U(0,1) | anchor (uniform, P-3) | IN |
| 1 | min(Z, 0) + 1% noise | **decision F1** | IN (R1 + R7) |
| 2 | min(Z, 0) + 5% noise | **decision F1** | IN (R1 + R7) |
| 3 | min(B, q) + 1% noise, a = 1.2 | **decision F2** | IN (R3 + R7) |
| 4 | min(B, q) + 5% noise, a = 1.2 | **decision F2** | IN (R3 + R7) |
| 5 | 0.1δ_0 + 0.9 Exp(1), + 1% noise | **decision F3** | IN (R4 + R7) |
| 6 | 0.1δ_0 + 0.9 Exp(1), + 5% noise | **decision F3** | IN (R4 + R7) |
| 7 | min(Z, 0) | atom row | IN (R1) |
| 8 | min(B, q), a = 1.2 | atom row | IN (R3) |
| 9 | 0.1δ_0 + 0.9 Exp(1) | atom row | IN (R4) |
| 10-12 | F1, F2, F3 + 10% noise | descriptive | IN (R7) |
| 13, 14 | F1 with the atom spread as U[0, δ], δ = 1%, 5% IQR_Z | descriptive (U-spread twin) | IN (R8: φ rising to 0, then plateau 0.5/δ > φ(0)) |
| 15, 16 | F2 with the atom spread as U[q, q + δ], δ = 1%, 5% IQR_B | descriptive | IN (R8: rising density, then plateau 0.1/δ ≥ density at q) |
| 17, 18 | F3 with the atom spread as U[0, δ], δ = 1%, 5% ln 3 | descriptive | IN (R8: 0.1/δ + 0.9e^{−x} on (0, δ), then 0.9e^{−x}; nonincreasing) |
| 19, 20 | min(Z, −0.5): atom; + 5% noise | descriptive | IN (R1, c ≤ 0; + R7) |
| 21, 22 | min(Z, −1): atom; + 5% noise | descriptive | IN |
| 23, 24 | min(U, 0.9) (a = 1): atom; + 5% noise (IQR 0.5) | descriptive (= atom version of H&H n = 4 witness family) | IN (R3) |
| 25, 26 | min(B, q), a = 1.5: atom; + 5% noise | descriptive | IN (R3) |
| 27, 28 | max(0, Z − 0.5): atom; + 5% noise | descriptive (ReLU, b = −0.5) | IN (R5) |
| 29 | W(0.1, 0.005) two-step | **known-excess anchor** (E-6) | IN (R8, M = [0.995, 1]) |
| 30-34 | W + slope β ∈ {0.1, 0.25, 0.5, 1, 2} (CV-9) | K10 | IN (R8: density nondecreasing; at L, shoulder ≤ 1.81 < 20 = top) |
| 35-37 | W with logistic-rounded top τ ∈ {0.003, 0.01, 0.03} (CV-9; shoulder mass 0.9 uniform on (0, 1)\*) | K10 | IN (sum of a constant and an increasing density) |
| 38 | unique-mode ramp, w = 0.1, ε_r = 0.01 (CV-9) | K10 / unique-mode witness | IN, M = {1} |
| 39 | off-edge pile-up, m0 = 0.8, w = 0.1, ε = 0.005 (CV-9) | K10 (= K2 interior) | IN, M = [0.795, 0.8] |
| 40-43 | W(w, 0.005), w ∈ {0.02, 0.05, 0.2, 0.4} | K10 | IN (top density w/ε ≥ 4 > 0.985) |
| 44 | t5 clipped to [−2, 2] | **non-null law** | OUT (R6) |
| 45 | max(0, Z + 0.5) | **non-null law** | OUT (R5, b > 0) |
| 46 | N(0,1) rounded to the 0.5 grid | **non-null law** | OUT (R6) |
| 47, 48 | W(0.1, 0.005) + σE, σ ∈ {0.005, 0.025} (1%, 5% of the shoulder IQR ≈ 0.5) | descriptive (noised witness) | IN (R8 + R7) |

\* f = 35-37: the shoulder component is U(0, 1) with mass 0.9, and the top component has mass 0.1 with density ∝ σ((x − 0.995)/τ) on (0, 1). The total density is nondecreasing on (0, 1).

## 4. Anchors (P-3, E-6)
- **f = 0 (uniform):** the CP interval for P(D ≥ q̂) must contain α at every (n, α), using the counts at q̂. A 3-SE band is the secondary check.
- **f = 29 (witness):** the point rejection at α = 0.05 must be ≥ 0.15 at n = 200. If either anchor fails, the screen is void and no decision is read.

## 5. Decision rule (adopted verbatim; α = 0.05; n ∈ {200, 1000}; families F1 = {1, 2}, F2 = {3, 4}, F3 = {5, 6})
- Simultaneous bounds over the 12 decision cells: one-sided CP at level 1 − 0.05/12 per cell, plus the bracket failure γ.
- **A family is MATERIAL** iff at least one of its 4 cells has a simultaneous lower bound (counts at X_(u)) ≥ 0.075.
- **A family is CONTROLLED** iff all 4 of its cells have a simultaneous upper bound (counts at X_(l)) < 0.075. Otherwise it is INCONCLUSIVE.
- **"Practically broad on the panel"** iff ≥ 2 families are MATERIAL.
- **Falsified on the panel** iff all 3 families are CONTROLLED.
- **Otherwise:** NOT ESTABLISHED, with the pattern reported.
- Pre-committed wording for the verdicts (a)-(d) was fixed before the run (not reproduced here).

## 6. What this screen does NOT establish
- A size (the families are named laws, and a search is not a supremum).
- Any asymptotic statement.
- Prevalence in real data.
- Anything from the descriptive rows beyond a description.

Outside-null rows are rejection under a non-null law (power), never size.

## 7. Recorded S predictions (before running)
- β ≤ 0.25 keeps a clear excess at n ≤ 1000; β ≥ 1 is near or below 5% by n = 500.
- F1 (censored at the peak, nearly flat just below it): a modest excess in the atom row; smaller with noise.
- F2 a = 1.2: an excess in the atom row; noise at 5% reduces it.
- F3: conservative (decreasing density, no flat part).
- The noised witness at 1%: the excess largely survives at n ≤ 200 and shrinks as n grows (σ fixed, so resolution grows). At 5%: much reduced.
- Decision: I expect "NOT ESTABLISHED" or "broad" to hinge on F1/F2 at 1% noise; F3 CONTROLLED.
