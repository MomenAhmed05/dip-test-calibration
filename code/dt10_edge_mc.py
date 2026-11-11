"""DT10 (Step 11, task 2): Brownian-limit check of C-edge via Lemma S/F (DT10-EDGE.md).

For a pair of independent discretised bridges b1, b2 (m cells; piecewise linear between grid values) and a
mode quantile u (a = sqrt(u), q = sqrt(1 - u)), Lemma F gives:
  d_u <= e  iff  e >= max(a G(b1), q G(b2))/2  and  a L_{b1}(e/a) + q L_{b2}(e/q) <= 0,
with G(b) = max(b - GCM b) and L_b(e) = max(-e, max_{s<x<1} [(b(x) - e) + (1 - x)(b(x) - b(s) - 2e)/(x - s)])
(sup over vertices, Lemma F(e)). Inner max over s: tangent query on the lower hull of the prefix (binary
search), so L costs O(m log m).
Discretisation: E_i = ||b_i - PL(b_i)|| is drawn exactly and independently (max of m scaled Kolmogorov
variables, as in qbar_cert.py). Since d is 1-Lipschitz, with E_u = max(a E1, q E2):
  feasible at t - E_u  => true d_u <= t (certain NON-exceedance);  infeasible at t + E_u => true d_u > t.
Common random numbers: the same (b1, b2, E1, E2) for every u on the grid (Corollary S'), u = 1 - eps^2.
Output per path and (t, u): code 0 = certain non-exceedance, 2 = certain exceedance, 1 = undecided,
plus the plain PL indicator (d_PL > t).
Usage: python dt10_edge_mc.py m R shard nshards root tag
"""
import json
import os
import sys
import time

import numpy as np
from numba import njit
from scipy.stats import kstwobign

EPS = (0.0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.45, 0.6, 0.7071067811865476)   # eps = sqrt(1-u)
TS = (0.73139, 0.85987)


@njit(cache=True)
def gcm_gap(s, b):
    m = len(s)
    hx = np.empty(m)
    hy = np.empty(m)
    k = 0
    for i in range(m):
        while k >= 2 and (hy[k - 1] - hy[k - 2]) * (s[i] - hx[k - 2]) >= (b[i] - hy[k - 2]) * (hx[k - 1] - hx[k - 2]):
            k -= 1
        hx[k] = s[i]
        hy[k] = b[i]
        k += 1
    g = 0.0
    j = 0
    for i in range(m):
        while j < k - 2 and hx[j + 1] <= s[i]:
            j += 1
        if hx[j + 1] == hx[j]:
            v = hy[j]
        else:
            v = hy[j] + (hy[j + 1] - hy[j]) * (s[i] - hx[j]) / (hx[j + 1] - hx[j])
        if b[i] - v > g:
            g = b[i] - v
    return g


@njit(cache=True)
def lval(s, b, e):
    """L_b(e) of Lemma F(c) for the piecewise-linear b with vertices (s, b), s[0] = 0, s[-1] = 1."""
    m = len(s)
    hx = np.empty(m)
    hy = np.empty(m)
    k = 0
    best = -e
    for j in range(m - 1):
        if k >= 1 and j >= 1:
            qx = s[j]
            qy = b[j] - 2.0 * e
            lo = 0
            hi = k - 1
            while lo < hi:        # first hull index whose slope to Q is >= that of the next one
                mid = (lo + hi) // 2
                s1 = (qy - hy[mid]) / (qx - hx[mid])
                s2 = (qy - hy[mid + 1]) / (qx - hx[mid + 1])
                if s1 >= s2:
                    hi = mid
                else:
                    lo = mid + 1
            sig = (qy - hy[lo]) / (qx - hx[lo])
            v = (b[j] - e) + (1.0 - s[j]) * sig
            if v > best:
                best = v
        while k >= 2 and (hy[k - 1] - hy[k - 2]) * (s[j] - hx[k - 2]) >= (b[j] - hy[k - 2]) * (hx[k - 1] - hx[k - 2]):
            k -= 1
        hx[k] = s[j]
        hy[k] = b[j]
        k += 1
    return best


@njit(cache=True)
def feasible(e, a, q, G1, G2, s, b1, b2):
    if e <= 0.0 or 2.0 * e < a * G1 or 2.0 * e < q * G2:
        return False
    if q == 0.0:
        return True
    return a * lval(s, b1, e / a) + q * lval(s, b2, e / q) <= 0.0


def bridge(rng, m):
    z = rng.standard_normal(m)
    w = np.concatenate([[0.0], np.cumsum(z)]) / np.sqrt(m)
    t = np.arange(m + 1) / m
    return w - t * w[-1]


def interp_err(rng, m):
    v = rng.uniform()
    return float(kstwobign.isf(-np.expm1(np.log(v) / m))) / np.sqrt(m)


def main():
    m, R, shard, ns = (int(x) for x in sys.argv[1:5])
    root = int(sys.argv[5])
    tag = sys.argv[6]
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "deposits", "dt10")
    try:
        os.makedirs(out, exist_ok=True)
    except OSError as err:
        sys.exit(f"cannot create {out}: {err}")
    rng = np.random.default_rng(np.random.SeedSequence([root, m, shard]))
    s = np.arange(m + 1) / m
    t0 = time.time()
    CODE, PL, GG = [], [], []
    for _ in range(shard, R, ns):
        b1 = bridge(rng, m)
        b2 = bridge(rng, m)
        E1 = interp_err(rng, m)
        E2 = interp_err(rng, m)
        G1, G2 = gcm_gap(s, b1), gcm_gap(s, b2)
        code, pl = [], []
        for t in TS:
            cr, pr = [], []
            for eps in EPS:
                a, q = np.sqrt(1.0 - eps * eps), eps
                Eu = max(a * E1, q * E2)
                if feasible(t - Eu, a, q, G1, G2, s, b1, b2):
                    c = 0
                elif not feasible(t + Eu, a, q, G1, G2, s, b1, b2):
                    c = 2
                else:
                    c = 1
                cr.append(c)
                pr.append(0 if feasible(t, a, q, G1, G2, s, b1, b2) else 1)
            code.append(cr)
            pl.append(pr)
        CODE.append(code)
        PL.append(pl)
        GG.append([G1, G2, E1, E2])
    res = {"m": m, "root": root, "shard": shard, "nshards": ns, "eps": EPS, "ts": TS, "code": CODE, "pl": PL,
           "G1G2E1E2": GG, "seconds": time.time() - t0}
    path = os.path.join(out, f"{tag}_m{m}_r{root}_s{shard}.json")
    try:
        with open(path + ".tmp", "w") as fh:
            json.dump(res, fh)
        os.replace(path + ".tmp", path)
    except OSError as err:
        sys.exit(f"cannot write {path}: {err}")
    print(f"m={m} shard={shard} paths={len(CODE)} {time.time() - t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
