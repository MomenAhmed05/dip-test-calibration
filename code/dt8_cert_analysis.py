"""DT8 certificate analysis (DT8-PREREG adopted). No random numbers are drawn here.

Part B (finite-n uniform validity table): for each n and alpha, fixed-sequence test over the cutoff grid
from the LARGEST cutoff down; a cutoff is certified iff every block count K_j <= b(N, alpha, d), where
b = max{x : P(Bin(N, alpha) <= x) <= d}, d = 0.01/48 (equivalent to the one-sided CP upper limit at
error d being <= alpha). Stop at the first failure; report the smallest cutoff certified.
Part A (qbar bracket): order statistics hi_(k), lo_(l) with P(Bin(R, 1-a) >= k) <= eps and
P(Bin(R, 1-a) <= l-1) <= eps, eps = 0.01/4; then the worst-case uniform-calibrated size
1 - F_M(2 qbar) bracketed with the certified F_M enclosure (fm_cert_h0.0002.npz).
Usage: python dt8_cert_analysis.py B root | A root tag
"""
import glob
import json
import os
import sys

import mpmath as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DEP = os.path.join(HERE, "..", "deposits")
NS = [50, 100, 200, 400, 800, 1600]
GRID = {0.05: [0.75, 0.74, 0.735, 0.7315], 0.01: [0.90, 0.88, 0.87, 0.86]}   # largest first
D_B = 0.01 / 48
EPS_A = 0.01 / 4
mp.mp.dps = 50


def binom_cdf(x, N, p):
    """P(Bin(N, p) <= x) via the regularised incomplete beta (exact to mp precision)."""
    if x < 0:
        return mp.mpf(0)
    if x >= N:
        return mp.mpf(1)
    p = mp.mpf(p)
    q = 1 - p
    term = q ** N                      # pmf(0); exact recurrence pmf(k+1) = pmf(k) (N-k)/(k+1) p/q
    total = term
    r = p / q
    for k in range(int(x)):
        term = term * (N - k) / (k + 1) * r
        total += term
    return total


def accept_count(N, alpha, d):
    """Largest x with P(Bin(N, alpha) <= x) <= d (-1 if none)."""
    lo, hi = -1, N
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if binom_cdf(mid, N, alpha) <= d:
            lo = mid
        else:
            hi = mid
    return lo


def part_b(root):
    out = {"root": root, "d": D_B, "rows": []}
    for n in NS:
        S, failed = [], 0
        for f in sorted(glob.glob(os.path.join(DEP, "dt8blk", f"cert_n{n}_r{root}_s*.json"))):
            r = _load(f)
            S += r["S"]
            failed += len(r.get("failed", []))
        if not S:
            continue
        S = np.array(S)
        N = len(S)
        for alpha, grid in GRID.items():
            b = accept_count(N, alpha, D_B)
            cert, trail = None, []
            for c in grid:
                K = [int(k) for k in np.greater(S, c).sum(axis=0)]
                ok = all(k <= b for k in K)
                trail.append({"c": c, "K": K, "pass": ok})
                if not ok:
                    break
                cert = c
            out["rows"].append({"n": n, "N": N, "alpha": alpha, "b": b, "failed_blocks": failed,
                                "certified_cutoff": cert, "trail": trail,
                                "edge_block_rate_at_cert": None if cert is None else
                                float((S[:, -1] > cert).mean())})
    _dump(out, os.path.join(DEP, "dt8blk", f"cert_analysis_r{root}.json"))
    for r in out["rows"]:
        print(f"n={r['n']:5d} N={r['N']:6d} alpha={r['alpha']}: b={r['b']}, certified c={r['certified_cutoff']}, "
              f"trail={[(t['c'], max(t['K']), t['pass']) for t in r['trail']]} failed={r['failed_blocks']}")


def part_a(root, tag):
    lo, hi = [], []
    for f in sorted(glob.glob(os.path.join(DEP, "qbar", f"{tag}_s*.json"))):
        r = _load(f)
        if r["root"] != root:
            raise ValueError(f"root mismatch in {f}")
        lo += r["lo"]
        hi += r["hi"]
    width = float(np.median(np.array(hi) - np.array(lo)))
    lo, hi = np.sort(lo), np.sort(hi)
    R = len(lo)
    fm = np.load(os.path.join(DEP, "dt8", "fm_cert_h0.0002.npz"))
    ys, flo, fhi = fm["ys"], fm["lo"], fm["hi"]
    res = {"root": root, "tag": tag, "R": R, "eps": EPS_A, "alphas": {}}
    for a in (0.05, 0.01):
        p = 1 - a
        k = _smallest(lambda kk, p=p: 1 - binom_cdf(kk - 1, R, p) <= EPS_A, R)    # P(Z >= k) <= eps
        ell = _largest(lambda ll, p=p: binom_cdf(ll - 1, R, p) <= EPS_A, R)       # P(Z <= l-1) <= eps
        if k is None or ell is None:
            print(f"alpha={a}: R={R} too small for eps={EPS_A}")
            continue
        q_up, q_lo = float(hi[k - 1]), float(lo[ell - 1])
        # size = 1 - F_M(2 qbar), F_M nondecreasing: bracket with the enclosure at outward grid points
        i_up = int(np.searchsorted(ys, 2 * q_up, side="left"))          # ys[i_up] >= 2 q_up
        i_lo = int(np.searchsorted(ys, 2 * q_lo, side="right")) - 1     # ys[i_lo] <= 2 q_lo
        size_lo = 1 - float(fhi[i_up])
        size_hi = 1 - float(flo[i_lo])
        res["alphas"][str(a)] = {"k": k, "l": ell, "qbar_bracket": [q_lo, q_up],
                                 "worst_size_bracket": [size_lo, size_hi],
                                 "median_width": width}
        print(f"alpha={a}: qbar in [{q_lo:.5f}, {q_up:.5f}] (l={ell}, k={k}, R={R}); "
              f"worst uniform-calibrated size in [{size_lo:.4f}, {size_hi:.4f}]")
    _dump(res, os.path.join(DEP, "qbar", f"{tag}_analysis.json"))


def _load(path):
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, ValueError) as exc:
        raise SystemExit(f"cannot read {path}: {exc}") from exc


def _dump(obj, path):
    try:
        with open(path, "w") as fh:
            json.dump(obj, fh, indent=1)
    except OSError as exc:
        raise SystemExit(f"cannot write {path}: {exc}") from exc


def _smallest(pred, R):
    """Smallest k in 1..R with pred(k), pred monotone False..True; None if none."""
    if not pred(R):
        return None
    lo, hi = 0, R
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if pred(mid):
            hi = mid
        else:
            lo = mid
    return hi


def _largest(pred, R):
    """Largest l in 1..R with pred(l), pred monotone True..False; None if none."""
    if not pred(1):
        return None
    lo, hi = 1, R + 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if pred(mid):
            lo = mid
        else:
            hi = mid
    return lo


if __name__ == "__main__":
    if sys.argv[1] == "B":
        part_b(int(sys.argv[2]))
    else:
        part_a(int(sys.argv[2]), sys.argv[3])
