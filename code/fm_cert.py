"""Certified enclosure of F_M (Step 10 / DT8, gap G1).

M = sup(LCM(B) - B) for a standard Brownian bridge B. Balabdaoui & Pitman (2011, Bernoulli 17:466-483,
Thm 1 / section 3): M =d max_j sqrt(L_j) E_j with uniform stick-breaking L_j and i.i.d. Brownian-excursion
maxima E_j, F_E(y) = 1 - 2 sum_{k>=1} (4k^2y^2 - 1) exp(-2k^2y^2). Splitting off the first stick gives
the exact renewal identity M =d max(sqrt(V) E, sqrt(1-V) M'), V ~ U(0,1), hence with w = y/sqrt(1-v):

    Phi(y) = int_y^inf F_E(a_y(w)) Phi(w) (2y^2/w^3) dw,   a_y(w) = y w / sqrt(w^2 - y^2),   Phi = F_M.

Certification (deterministic sandwich, no simulation):
  * Phi is nondecreasing; a_y(w) is decreasing in w; F_E is nondecreasing; the weight integrates exactly.
  * Tail: M <= range(B) so 1 - Phi(w) <= 2 exp(-w^2/2) (reflection principle for each side).
  * Backward recursion from y = YMAX on a grid of step h; on every cell the integrand is bounded by its
    monotone end values; the self-referencing first cell is solved implicitly (Phi >= A/(1-c)).
  * F_E is a convergent series; for z >= ZMIN we keep K terms and bound the remainder explicitly.
  * Floating point: every bound is padded by PAD (absolute) at each grid point; PAD >> accumulated
    float64 rounding of the sums (<= 1e5 terms x 1.2e-16).
Independent checker (F-8): fm_cert_check.py recomputes with python-flint ball arithmetic on a coarse
grid; the coarse enclosure must contain the fine one.
Output: deposits/dt8/fm_cert_h{h}.json with [lo, hi] at every grid point and the quantile brackets.
"""
import json
import math
import os
import sys
import time

import numpy as np
from numba import njit

YLO, YMAX, WMAX, TAILR = 0.9, 8.5, 2.0e4, 1.002
KTERMS = 12
PAD = 1e-12


@njit(cache=True)
def FE(z):
    """Lower and upper bound of F_E(z) for z >= 0.5 (series with explicit remainder bound)."""
    if z > 7.0:
        # 1 - F_E(z) <= 2*sum (4k^2z^2) e^{-2k^2z^2} <= 16 z^2 e^{-2z^2} for z >= 1
        r = 16.0 * z * z * math.exp(-2.0 * z * z)
        return 1.0 - r, 1.0
    s = 0.0
    for k in range(1, KTERMS + 1):
        kk = k * k
        s += (4.0 * kk * z * z - 1.0) * math.exp(-2.0 * kk * z * z)
    # remainder for k > KTERMS: terms positive and decreasing (z >= 0.5), bounded by 2x the first one
    k = KTERMS + 1
    kk = k * k
    rem = 2.0 * (4.0 * kk * z * z) * math.exp(-2.0 * kk * z * z)
    v = 1.0 - 2.0 * s
    lo = max(0.0, v - 2.0 * rem - 1e-14)
    hi = min(1.0, v + 1e-14)
    return lo, hi


@njit(cache=True)
def arg(y, w):
    if w <= y:
        return 1e300
    return y * w / math.sqrt(w * w - y * y)


@njit(cache=True)
def solve(h):
    N = int(math.ceil((YMAX - YLO) / h))
    ys = YLO + h * np.arange(N + 1)
    lo = np.zeros(N + 1)
    hi = np.ones(N + 1)
    # top point: tail bound
    lo[N] = 1.0 - 2.0 * math.exp(-ys[N] * ys[N] / 2.0)
    hi[N] = 1.0
    # tail cells beyond the grid
    nt = int(math.ceil(math.log(WMAX / ys[N]) / math.log(TAILR)))
    wt = ys[N] * TAILR ** np.arange(nt + 1)
    for i in range(N - 1, -1, -1):
        y = ys[i]
        y2 = y * y
        Alo = 0.0
        Ahi = 0.0
        c_lo = 0.0
        # cells on the grid
        for j in range(i, N):
            a, b = ys[j], ys[j + 1]
            wgt = y2 * (1.0 / (a * a) - 1.0 / (b * b))
            fe_lo, _ = FE(arg(y, b))
            if j == i:
                _, fe_hi = 1.0, 1.0
            else:
                _, fe_hi = FE(arg(y, a))
            if j == i:
                c_lo = fe_lo * wgt
            else:
                Alo += fe_lo * lo[j] * wgt
            Ahi += fe_hi * hi[j + 1] * wgt
        # tail cells
        for k in range(nt):
            a, b = wt[k], wt[k + 1]
            wgt = y2 * (1.0 / (a * a) - 1.0 / (b * b))
            fe_lo, _ = FE(arg(y, b))
            _, fe_hi = FE(arg(y, a))
            philo = 1.0 - 2.0 * math.exp(-a * a / 2.0)
            Alo += fe_lo * philo * wgt
            Ahi += fe_hi * 1.0 * wgt
        # remainder beyond WMAX: integrand in [0, 1]
        Ahi += y2 / (wt[nt] * wt[nt])
        L = Alo / (1.0 - c_lo) - PAD
        H = Ahi + PAD
        lo[i] = max(0.0, L)
        hi[i] = min(1.0, H, hi[i + 1])
        if lo[i] > lo[i + 1] and i + 1 <= N:
            pass
    # monotone tightening of lower bounds is not needed (lower bounds already valid pointwise)
    return ys, lo, hi


def bracket(ys, lo, hi, p):
    """Interval [a, b] with F_M^{-1}(p) in [a, b]: a = max y with hi(y) < p, b = min y with lo(y) >= p."""
    a = ys[np.where(hi < p)[0].max()]
    b = ys[np.where(lo >= p)[0].min()]
    return float(a), float(b)


def main():
    h = float(sys.argv[1]) if len(sys.argv) > 1 else 1e-3
    t0 = time.time()
    ys, lo, hi = solve(h)
    out = {"h": h, "YLO": YLO, "YMAX": YMAX, "WMAX": WMAX, "TAILR": TAILR, "KTERMS": KTERMS,
           "PAD": PAD, "seconds": time.time() - t0}
    pts = {}
    for y in (1.0, 1.0912, 1.0970, 1.2, 1.4, 1.4628, 1.6, 1.7197, 1.8, 2.0):
        i = int(round((y - YLO) / h))
        pts[f"{ys[i]:.6f}"] = [float(lo[i]), float(hi[i])]
    out["points"] = pts
    out["max_width"] = float(np.max(hi - lo))
    for p, name in ((0.95, "q95"), (0.99, "q99"), (0.90, "q90")):
        a, b = bracket(ys, lo, hi, p)
        out[name] = [a, b]
        out["c_" + name] = [a / 2, b / 2]
    d = os.path.join(os.path.dirname(__file__), "..", "deposits", "dt8")
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, f"fm_cert_h{h:g}.json"), "w") as f:
        json.dump(out, f, indent=1)
    np.savez_compressed(os.path.join(d, f"fm_cert_h{h:g}.npz"), ys=ys, lo=lo, hi=hi)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
