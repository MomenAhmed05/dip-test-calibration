"""DT10 DEV checks (root 2026110730; not CONFIRM).
(1) lval (hull + binary search) == brute-force O(m^2) formula of Lemma F(c), random PL bridges.
(2) Lemma F as a whole: d_u by bisection on `feasible` == direct LP distance of the PL path H_u to K_1/2.
(3) timing of one path (all u, both t) for m = 2^14 and 2^16.
"""
import json
import os
import sys
import time

import numpy as np
from scipy.optimize import linprog

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dt10_edge_mc as E  # noqa: E402

ROOT = 2026110730


def l_brute(s, b, e):
    best = -e
    for j in range(1, len(s) - 1):
        for i in range(j):
            v = (b[j] - e) + (1 - s[j]) * (b[j] - b[i] - 2 * e) / (s[j] - s[i])
            best = max(best, v)
    return best


def d_lp(b1, b2, u):
    """Direct LP: distance of the PL path H_u (Lemma S) to K_1/2, on the 2m+1 knots."""
    m = len(b1) - 1
    a, q = np.sqrt(u), np.sqrt(1 - u)
    x = np.concatenate([np.arange(m + 1) / (2 * m), 0.5 + np.arange(1, m + 1) / (2 * m)])
    r2 = -b2[::-1]
    h = np.concatenate([a * b1, q * r2[1:]])
    K = len(x)
    A, rhs = [], []
    for j in range(K):
        r = np.zeros(K + 1)
        r[j], r[K] = 1, -1
        A.append(r)
        rhs.append(h[j])
        r = np.zeros(K + 1)
        r[j], r[K] = -1, -1
        A.append(r)
        rhs.append(-h[j])
    for j in range(1, K - 1):
        h1, h2 = x[j] - x[j - 1], x[j + 1] - x[j]
        r = np.zeros(K + 1)
        if j < m:          # convex on [0, 1/2] (interior knots)
            r[j - 1], r[j], r[j + 1] = -1 / h1, 1 / h1 + 1 / h2, -1 / h2
        elif j > m:        # concave on [1/2, 1]
            r[j - 1], r[j], r[j + 1] = 1 / h1, -1 / h1 - 1 / h2, 1 / h2
        else:
            continue
        A.append(r)
        rhs.append(0.0)
    c = np.zeros(K + 1)
    c[K] = 1
    res = linprog(c, A_ub=np.array(A), b_ub=np.array(rhs), bounds=[(None, None)] * K + [(0, None)],
                  method="highs")
    return float(res.x[K])


def d_bisect(s, b1, b2, u):
    a, q = np.sqrt(u), np.sqrt(1 - u)
    G1, G2 = E.gcm_gap(s, b1), E.gcm_gap(s, b2)
    lo, hi = 0.0, 5.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if E.feasible(mid, a, q, G1, G2, s, b1, b2):
            hi = mid
        else:
            lo = mid
    return hi


def main():
    rng = np.random.default_rng(np.random.SeedSequence([ROOT, 1]))
    out = {}
    worst = 0.0
    for _ in range(30):
        m = int(rng.integers(5, 120))
        s = np.arange(m + 1) / m
        b = E.bridge(rng, m)
        G = E.gcm_gap(s, b)
        for e in (G / 2, G / 2 + 0.01, G, 3 * G + 0.5):
            worst = max(worst, abs(E.lval(s, b, e) - l_brute(s, b, e)))
    out["lval_vs_brute_maxabs"] = worst
    worst, cases = 0.0, []
    for _ in range(40):
        m = int(rng.integers(4, 60))
        s = np.arange(m + 1) / m
        b1, b2 = E.bridge(rng, m), E.bridge(rng, m)
        for u in (0.5, 0.7, 0.9, 0.99, 1.0 - 1e-6, 1.0):
            dl, db = d_lp(b1, b2, u), d_bisect(s, b1, b2, u)
            worst = max(worst, abs(dl - db))
            if abs(dl - db) > 1e-7:
                cases.append([m, u, dl, db])
    out["lemmaF_vs_LP_maxabs"] = worst
    out["lemmaF_vs_LP_bad"] = cases[:10]
    tm = {}
    for m in (2 ** 14, 2 ** 16):
        s = np.arange(m + 1) / m
        b1, b2 = E.bridge(rng, m), E.bridge(rng, m)
        G1, G2 = E.gcm_gap(s, b1), E.gcm_gap(s, b2)
        t0 = time.time()
        for t in E.TS:
            for eps in E.EPS:
                a, q = np.sqrt(1 - eps * eps), eps
                for ee in (t - 0.01, t, t + 0.01):
                    E.feasible(ee, a, q, G1, G2, s, b1, b2)
        tm[m] = time.time() - t0
    out["sec_per_path"] = tm
    p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "deposits", "dt10", "dev_check.json")
    try:
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, "w") as fh:
            json.dump(out, fh, indent=1)
    except OSError as err:
        sys.exit(f"cannot write {p}: {err}")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
