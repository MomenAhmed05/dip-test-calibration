"""Step 6 (DT4) independent check of F_M, the CDF of M = sup(LCM(B) - B) for a Brownian bridge B.

Method (independent of the Balabdaoui-Filali numerical inversion used in dt1b_confirm.bf_table):
Balabdaoui & Pitman (2011, Bernoulli 17(1), Thm 1 and eq. (5)): M =d max_j sqrt(L_j) M3_j, with (L_j) uniform
stick-breaking and M3_j i.i.d. maxima of a standard Brownian excursion, so F_M(x) = E prod_j F3(x / sqrt(L_j)),
F3(y) = 1 + 2 sum_{k>=1} (1 - 4k^2 y^2) exp(-2k^2 y^2)  (excursion maximum = Kuiper law).
The expectation is over the stick-breaking only (the integrand is exact), so the only random error is Monte Carlo
over L, bounded by Hoeffding (integrand in [0, 1]). Deterministic errors (series truncation at K terms, stick
truncation after J pieces) are bounded below and added.

Root 2026110640 (Step 6 numerics; registered in CONVENTIONS). Run: python code/fm_check.py [--quick]
"""

import json
import math
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dt1b_confirm import bf_table  # noqa: E402

QUICK = "--quick" in sys.argv
ROOT = 2026110640
N = 200_000 if QUICK else 10_000_000
CHUNK = 100_000
J = 60  # stick pieces kept; remainder bounded below
K = 10  # series terms (tail <= 1e-60 for y >= 0.6)
YCUT = 6.0  # pieces with x/sqrt(L) >= YCUT get factor 1; 1 - F3(y) <= 8 y^2 e^{-2y^2} <= 1.5e-29 there
XS = np.round(np.arange(0.60, 2.001, 0.02), 4)
DELTA = 1e-6  # Hoeffding failure probability per point
OUT = Path(__file__).resolve().parents[1] / "deposits" / ("dt4_quick" if QUICK else "dt4")


def f3(y):
    """Excursion-maximum CDF for y >= 0.5 (series converges fast there); y may be large."""
    y = np.minimum(y, 50.0)
    s = np.zeros_like(y)
    for k in range(1, K + 1):
        a = 2.0 * k * k * y * y
        s += (1.0 - 2.0 * a) * np.exp(-a)
    return 1.0 + 2.0 * s


def series_tail_bound(ymin):
    """|sum_{k>K} 2(1 - 4k^2 y^2) e^{-2k^2 y^2}| <= sum_{k>K} 8 k^2 y^2 e^{-2k^2 y^2} for y >= ymin >= 0.5."""
    t = 0.0
    for k in range(K + 1, K + 200):
        a = 2.0 * k * k * ymin * ymin
        t += 4.0 * a * math.exp(-a)
    return t


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT)))
    sums = np.zeros(len(XS))
    max_rem = 0.0
    for _ in range(N // CHUNK):
        w = rng.random((CHUNK, J))
        rem = np.cumprod(1.0 - w, axis=1)
        L = w * np.concatenate([np.ones((CHUNK, 1)), rem[:, :-1]], axis=1)
        max_rem = max(max_rem, float(rem[:, -1].max()))
        sq = np.sqrt(L)
        for i, x in enumerate(XS):
            y = x / sq
            m = y < YCUT
            fac = np.ones_like(y)
            fac[m] = f3(y[m])
            sums[i] += np.prod(fac, axis=1).sum()
    est = sums / N
    hoeff = math.sqrt(math.log(2 / DELTA) / (2 * N))
    # remainder pieces (total length R <= max_rem): sum of (1 - F3) over them <= 8 x^2 e^{-2x^2/R} / R (see DT4-CERT)
    rem_bound = max(8 * x * x * math.exp(-2 * x * x / max_rem) / max_rem for x in XS) if max_rem > 0 else 0.0
    # |prod a_j - prod b_j| <= sum |a_j - b_j| for factors in [0, 1]: J pieces, each with series error <= tail
    ser = J * series_tail_bound(float(XS.min())) + J * 8 * YCUT**2 * math.exp(-2 * YCUT**2)
    xs_bf, fv_bf = bf_table()
    rows = []
    for x, e in zip(XS, est, strict=True):
        bf = float(np.interp(x, xs_bf, fv_bf)) if xs_bf[0] <= x <= xs_bf[-1] else None
        rows.append({"x": float(x), "F_M_est": float(e), "lo": float(e - hoeff - ser - rem_bound),
                     "hi": float(e + hoeff + ser), "BF_table": bf,
                     "diff_vs_BF": None if bf is None else float(e - bf)})
    out = {"root": ROOT, "N": N, "J": J, "K": K, "delta_per_point": DELTA, "hoeffding_eps": hoeff,
           "series_tail_bound": ser, "max_remainder_stick": max_rem, "remainder_bound": rem_bound,
           "rows": rows, "seconds": time.time() - t0}
    (OUT / "fm_check.json").write_text(json.dumps(out, indent=1))
    print(f"N={N} hoeff={hoeff:.2e} series={ser:.1e} max_rem={max_rem:.1e} rem_bound={rem_bound:.1e}")
    show = {0.6, 0.8, 1.0, 1.08, 1.1, 1.14, 1.16, 1.2, 1.26, 1.34, 1.4, 1.46, 1.6, 1.8, 2.0}
    for r in rows:
        if round(r["x"], 2) in show:
            print(f"x={r['x']:.2f} F_M={r['F_M_est']:.5f} [{r['lo']:.5f},{r['hi']:.5f}] BF={r['BF_table']} "
                  f"diff={r['diff_vs_BF']}")
    print(f"done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
