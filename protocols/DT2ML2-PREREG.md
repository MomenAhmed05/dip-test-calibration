# DT2ML2-PREREG — Step 4b part 2: longer, finer adversary training. DEV search; fresh re-estimates.

Frozen 2026-11-05 (file last modified 19:50:00 UTC), before any part-2 output (result file last modified 19:52:40 UTC). Reason: part 1 stopped before converging (still rising). Its 64-bin grid cannot represent pile-ups narrower than 1/64, while the best 2-block law has ε ≈ 0.011.
Roots (CV-6): **2026110623** training, **2026110624** re-estimation.

## Changes from part 1 (everything else as `DT2ML-PREREG.md`)
- K = 256 bins.
- 3,000 Adam steps, lr 0.03, batch 2,048 at n = 200.
- Smoothing τ annealed geometrically from 0.2 to 0.02.
- Runs (all unrestricted):
  - peak bin at 255 (right edge), 224, 192 (x = 0.75), 160, 128 (centre), from a uniform start;
  - **centre-asymmetric start** (peak bin 128; left steps ≈ 0, right steps steep): tests if the excess is reachable with a central peak;
  - one curved run (δ = 0.05 per bin in log-density... scaled to K = 256: δ = 0.0125 per bin, the same total slope as part 1), peak bin 255, uniform start.
- Re-estimation: 100,000 fresh datasets at n = 200 and 1000 for each final law.
- Baseline: the part-1 2-block optimum (fresh 25.75% at n = 200, 24.66% at n = 1000).

## Readings fixed in advance
- **Q-b2:** a learned law beats the 2-block baseline at n = 200 iff its fresh rejection minus 25.75% is > 2 × its binomial SE (≈ 0.28 points) → "ADVERSARY BEATS 2-BLOCK". Otherwise "2-BLOCK NOT BEATEN".
- **Q-c2:** the centre-asymmetric run is > 7.5% at n = 200 → "CENTRAL PEAK CAN FAIL" (consistent with BF: one heavy shoulder).
- **Q-d2:** the curved run ≤ 5.5% → "FLATS NEEDED" confirmed.

## Recorded S predictions
- The right-edge run reaches 25-26% with a shape close to the 2-block optimum (shoulder mass ≈ 0.89, pile-up ≈ 0.11 in about 3 bins).
- Q-b2: 2-BLOCK NOT BEATEN (within 1 point).
- Q-c2: yes, about 20%.
- Q-d2: FLATS NEEDED.
