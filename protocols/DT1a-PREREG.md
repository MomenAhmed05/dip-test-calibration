# DT1a-PREREG — Step 3a mechanism + convention audit. DEV-mechanism tier.

Frozen 2026-11-05 (file last modified 18:54:45 UTC), before any 3a output (result file last modified 18:55:28 UTC).

**Tier:** DEV-mechanism. This stage decides which formula the confirmation (3b) tests. It confirms nothing.

## 1. Laws (CV-9)
- f0: W(0.1, 0.005). Shoulder (0, 0.995) with mass 0.9; flat top (0.995, 1).
- f1: W(0.4, 0.005). Shoulder mass 0.6.
- f2: the unique-mode ramp (w = 0.1, ε_r = 0.01). Shoulder (0, 0.99) with mass 0.9, then a strictly rising density on (0.99, 1); M = {1}.

n ∈ {800, 3200, 12800}; **10,000 datasets per cell**. Root **2026110560** (DEV, declared in CV-6), spawn (f, n).

## 2. Per-sample quantities
- **D** = D_alg (CV-1).
- **V_sh** = ½·max over the shoulder points (x < L, where L = 0.995 for W and 0.99 for the ramp) of [F_n(x) − GCM of the lower corners] = the distance of F_n on the shoulder to convex functions (`path_smoke.vn_sample`).
- **T_top** = (k_top/n)·D_alg(the top points only), where k_top = the number of points ≥ L. This is S's proxy for the plateau's own free-mode fit cost. For the ramp it is computed the same way.
- The signed **√n(D − V_sh)**; equality |D − V_sh| ≤ 1e-12; the event D > V_sh + 1e-12.
- **√n·max(V_sh, T_top)** vs √n·D.

## 3. Tests against the Balabdaoui-Pitman law (bypass q̄)
- Z_n = 2√n·D/√(1−w) and Y_n = 2√n·V_sh/√(1−w). Compare their empirical CDFs with F_M at x ∈ {0.8, 1.0, 1.2, 1.4, 1.6}.
- F_M (Balabdaoui & Filali 2011, Table 2): 0.220644, 0.560962, 0.806055, 0.928968, 0.978091.
- **CONSISTENT at n** iff all 5 |ECDF − F_M| ≤ 0.015. (The MC SE is ≤ 0.005; the allowance 0.015 covers the finite-n bias seen in SM-6.)
- **Reading per law:**
  - "formula holds" iff Z_n is CONSISTENT at n = 12800, and its maximum discrepancy at 12800 is ≤ that at 800 + 0.005;
  - "shoulder-only fails" iff Y_n is CONSISTENT at 12800 but Z_n is not;
  - "both fail": then audit V_sh and the boundary conventions.

## 4. Convention fixtures (deterministic)
- all-identical (n = 5) → D_alg = 0 (CV-1);
- n = 2 and n = 3 distinct → 1/(2n) (R documentation: "for n ≤ 3, D_n is always 1/(2n)");
- ties on both jump sides (2-2-1 pattern), a mode atom (uniform + atom at the max);
- reflection x → −x and an affine map 3x + 7: equal D to 1e-15;
- strict vs non-strict rejection: the number of reference dips equal to q̂ (DT0-d references) = 0, so D > c and D ≥ c agree.

## 5. Recorded S predictions (before running)
- **f0 W(0.1):** Z_n CONSISTENT at all n; P(D = V_sh) ≥ 0.95 at n = 12800; signed √n(D − V_sh) has q99 ≤ 0.01. **Formula holds.**
- **f1 W(0.4):** Y_n CONSISTENT, but **Z_n not CONSISTENT**: D exceeds V_sh noticeably, because √0.4·U (plateau) is comparable to √0.6·M/2. The max(V_sh, T_top) tracks D much better than V_sh. **Shoulder-only fails at w = 0.4.**
- **f2 ramp:** like f0. Formula holds; P(D = V_sh) slightly lower than f0.
- **Fixtures:** all pass.

## 6. What this stage does NOT establish
- Any confirmation.
- The limit (n ≤ 12800).
- q̄_α.
- The proof of the reduction (Step 5).
