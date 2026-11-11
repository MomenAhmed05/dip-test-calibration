"""Certified bracket for qbar_alpha = (1 - alpha) quantile of D(B) (Step 10 / DT8, gap G1).

Per Brownian-bridge path (exact Gaussian grid with m cells):
  1. B_m = linear interpolation of the bridge at j/m.
  2. D(B_m) exactly up to lam/(2N): g -> a*g + l (a > 0, l affine on [0,1], constant outside) maps the H&H
     class U* onto itself (H&H 1985 Thm 2; U* is NOT a convex cone), so D(a f + l) = a D(f) and
     D(B_m) = lam * dip(H), H = id + B_m/lam, a strictly increasing
     piecewise-linear CDF on [0, 1] when lam > max |slope of B_m|. The ECDF of the N points
     H^{-1}((i - 1/2)/N) is within 1/(2N) of H in sup norm, and the dip is 1-Lipschitz (H&H Thm 1 for the
     identification with rho(., U*)), so |lam * dip(ECDF) - D(B_m)| <= lam/(2N). dip(ECDF) by diptest.
  3. ||B - B_m|| = max_j sup|residual bridge on cell j| / sqrt(m); the residuals are iid standard bridges
     independent of the grid, so E_m is drawn from P(E <= x) = K(x sqrt(m))^m (K = law of sup|bridge|).
  4. Coupling: the skeleton and E_m admit a coupling with a standard Brownian bridge B under
     which lo = max(0, lam(dip - eta) - E_m) <= D(B) <= hi = lam(dip + eta) + E_m almost surely,
     eta = 1/(2N) + 1e-12 (float slack). Floating point is NOT interval-checked (stated caveat).
Output: deposits/qbar/{tag}_s{shard}.json with lo/hi/raw. Quantile bracket by order statistics (separate step).
Usage: python qbar_cert.py m N R shard nshards root tag
"""
import json
import os
import sys
import time

import diptest
import numpy as np
from scipy.stats import kstwobign


def one_path(rng, m, N):
    z = rng.standard_normal(m)
    w = np.concatenate([[0.0], np.cumsum(z)]) / np.sqrt(m)
    t = np.arange(m + 1) / m
    b = w - t * w[-1]
    slopes = np.diff(b) * m
    lam = 1.05 * float(np.max(np.abs(slopes))) + 1e-6   # > max |slope|, so H is strictly increasing
    h = t + b / lam
    levels = (np.arange(1, N + 1) - 0.5) / N
    x = np.interp(levels, h, t)
    res = diptest.dipstat(x)
    d = float(res[0] if isinstance(res, tuple) else res)
    v = rng.uniform()
    p_tail = -np.expm1(np.log(v) / m)          # P(single cell sup > x) such that max law is exact
    e = float(kstwobign.isf(p_tail)) / np.sqrt(m)
    raw = lam * d
    eta = 1.0 / (2 * N) + 1e-12
    return max(0.0, raw - lam * eta - e), raw + lam * eta + e, raw, e, lam * eta


def main():
    m, N, R, shard, ns, root = (int(a) for a in sys.argv[1:7])
    tag = sys.argv[7]
    out = os.path.join(os.path.dirname(__file__), "..", "deposits", "qbar")
    os.makedirs(out, exist_ok=True)
    rng = np.random.default_rng(np.random.SeedSequence([root, m, N, shard]))
    rows = []
    t0 = time.time()
    for _ in range(shard, R, ns):
        rows.append(one_path(rng, m, N))
    a = np.array(rows)
    res = {"m": m, "N": N, "root": root, "shard": shard, "nshards": ns, "lo": a[:, 0].tolist(),
           "hi": a[:, 1].tolist(), "raw": a[:, 2].tolist(), "E": a[:, 3].tolist(),
           "ecdf_err": a[:, 4].tolist(), "seconds": time.time() - t0}
    with open(os.path.join(out, f"{tag}_s{shard}.json"), "w") as f:
        json.dump(res, f)
    print(f"{tag} shard={shard} paths={len(rows)} {time.time() - t0:.1f}s "
          f"mean width {np.mean(a[:, 1] - a[:, 0]):.5f}")


if __name__ == "__main__":
    main()
