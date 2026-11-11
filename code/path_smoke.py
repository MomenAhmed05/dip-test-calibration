"""Step 2 path-functional normalisation smoke (PLAN Step 2 item 3, G-003). DEV tier.

V(f) := d_sup(f, convex functions) on an interval. Identity checked: V(f) = (1/2) max(f - GCM(f)).
(A) identity on random grid functions: LP value vs (1/2)·hull gap.
(B) V(B) for Brownian bridge B on grids m = 1000 / 4000 / 16000: quantiles and grid bias.
(C) sample inequality: D_alg(X) <= V_n(X) := (1/2) max_i (i/n - GCM_lower(x_i)) (distance of F_n to convex,
    i.e. mode-at-right, CDFs), checked on uniform samples.
(D) witness reduction: for W(0.1, 0.005) samples, sqrt(n)·D_n vs sqrt(n)·V_n restricted to the shoulder,
    and the predicted limit sqrt(1-w)·V(B) against the DT0-d witness dips at n = 12800.
Seeds: root 2026110542 (CV-6, DEV).
"""

import json
import sys
import time
from pathlib import Path

import numpy as np
from numba import njit
from scipy.optimize import linprog

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dt0c_screen import dips_of, g_uniform, g_witness  # noqa: E402, I001

ROOT = 2026110542
W = Path(__file__).resolve().parents[1]
OUT = W / "deposits" / "dt0b"
QS = (0.5, 0.9, 0.95, 0.99)


@njit(cache=True)
def hull_gap(x, y):
    """max_i (y_i - H(x_i)), H = greatest convex minorant of the points (x_i, y_i), x increasing."""
    m = x.shape[0]
    hx = np.empty(m)
    hy = np.empty(m)
    k = 0
    for i in range(m):
        while k >= 2 and (hx[k - 1] - hx[k - 2]) * (y[i] - hy[k - 2]) - (hy[k - 1] - hy[k - 2]) * (x[i] - hx[k - 2]) <= 0:
            k -= 1
        hx[k] = x[i]
        hy[k] = y[i]
        k += 1
    best = 0.0
    j = 0
    for i in range(m):
        while j < k - 2 and hx[j + 1] < x[i]:
            j += 1
        if hx[j + 1] == hx[j]:
            h = hy[j]
        else:
            h = hy[j] + (hy[j + 1] - hy[j]) * (x[i] - hx[j]) / (hx[j + 1] - hx[j])
        if y[i] - h > best:
            best = y[i] - h
    return best


@njit(cache=True)
def vn_sample(xs, n_total, offset):
    """(1/2) max over upper corners of (F_n - GCM of lower corners), for sorted points xs whose ranks start
    at offset + 1 in a sample of size n_total."""
    m = xs.shape[0]
    lo = np.empty(m)
    for i in range(m):
        lo[i] = (offset + i) / n_total
    # GCM of lower corners, evaluated at the same x; the upper corner adds 1/n_total
    gap = hull_gap(xs, lo)  # gap of lower corners themselves (normally 0)
    # recompute hull and the upper-corner gap
    hx = np.empty(m)
    hy = np.empty(m)
    k = 0
    for i in range(m):
        while k >= 2 and (hx[k - 1] - hx[k - 2]) * (lo[i] - hy[k - 2]) - (hy[k - 1] - hy[k - 2]) * (xs[i] - hx[k - 2]) <= 0:
            k -= 1
        hx[k] = xs[i]
        hy[k] = lo[i]
        k += 1
    best = 0.0
    j = 0
    for i in range(m):
        while j < k - 2 and hx[j + 1] < xs[i]:
            j += 1
        if hx[j + 1] == hx[j]:
            h = hy[j]
        else:
            h = hy[j] + (hy[j + 1] - hy[j]) * (xs[i] - hx[j]) / (hx[j + 1] - hx[j])
        g = lo[i] + 1.0 / n_total - h
        if g > best:
            best = g
    return 0.5 * max(best, gap)


def lp_convex_distance(f):
    m = f.shape[0]
    # variables g_0..g_{m-1}, t ; minimise t
    c = np.zeros(m + 1)
    c[-1] = 1.0
    rows, rhs = [], []
    for i in range(m):  # f - g <= t ; g - f <= t
        r = np.zeros(m + 1)
        r[i], r[-1] = -1.0, -1.0
        rows.append(r)
        rhs.append(-f[i])
        r = np.zeros(m + 1)
        r[i], r[-1] = 1.0, -1.0
        rows.append(r)
        rhs.append(f[i])
    for i in range(1, m - 1):  # convexity on a uniform grid: -(g_{i-1} - 2 g_i + g_{i+1}) <= 0
        r = np.zeros(m + 1)
        r[i - 1], r[i], r[i + 1] = -1.0, 2.0, -1.0
        rows.append(r)
        rhs.append(0.0)
    res = linprog(c, A_ub=np.array(rows), b_ub=np.array(rhs), bounds=[(None, None)] * (m + 1), method="highs")
    return res.fun


def main():
    t0 = time.time()
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(ROOT)))
    out = {"tier": "DEV", "root": ROOT}

    # (A) identity
    diffs = []
    grid = np.linspace(0, 1, 40)
    for _ in range(200):
        f = np.cumsum(rng.standard_normal(40))
        diffs.append(abs(lp_convex_distance(f) - 0.5 * hull_gap(grid, f)))
    out["A_identity"] = {"cases": 200, "max_abs_diff": float(max(diffs))}

    # (B) V(B) for Brownian bridges
    out["B_bridge"] = {}
    vref = None
    for m in (1000, 4000, 16000):
        N = 20000
        t = np.arange(m + 1) / m
        v = np.empty(N)
        for b in range(0, N, 1000):
            inc = rng.standard_normal((1000, m)) / np.sqrt(m)
            w = np.concatenate([np.zeros((1000, 1)), np.cumsum(inc, axis=1)], axis=1)
            br = w - t[None, :] * w[:, -1:]
            for i in range(1000):
                v[b + i] = 0.5 * hull_gap(t, br[i])
        out["B_bridge"][m] = {"N": N, "mean": float(v.mean()),
                              "quantiles": [float(x) for x in np.quantile(v, QS)]}
        vref = v
    # (C) D_alg <= V_n on uniform samples
    worst = -1.0
    ratio = []
    for _ in range(5):
        x = np.sort(g_uniform(rng, (1000, 2000)), axis=1)
        d = dips_of(x)
        for i in range(1000):
            vn = vn_sample(x[i], 2000, 0)
            worst = max(worst, d[i] - vn)
            ratio.append(d[i] / vn)
    out["C_dip_le_Vn_uniform"] = {"samples": 5000, "n": 2000, "max(D - V_n)": float(worst),
                                  "holds": bool(worst <= 1e-12),
                                  "median D/V_n": float(np.median(ratio))}

    # (D) witness reduction
    n = 12800
    gen = g_witness()
    x = np.sort(gen(rng, (2000, n)), axis=1)
    d = dips_of(x)
    L = 0.995
    vsh = np.empty(2000)
    for i in range(2000):
        k = int(np.searchsorted(x[i], L))
        vsh[i] = vn_sample(x[i, :k].copy(), n, 0)
    sd, sv = np.sqrt(n) * d, np.sqrt(n) * vsh
    dt0d = np.sqrt(n) * np.load(W / "deposits" / "dt0d" / "dips_f0_n12800.npy")
    ref = np.sqrt(n) * np.load(W / "deposits" / "dt0d" / "ref_n12800.npy")
    q95u = float(np.quantile(ref, 0.95, method="inverted_cdf"))
    pred = np.sqrt(0.9) * vref
    out["D_witness"] = {
        "n": n, "samples": 2000,
        "corr(sqrt n D, sqrt n V_sh)": float(np.corrcoef(sd, sv)[0, 1]),
        "median D/V_sh": float(np.median(d / vsh)), "frac D<=V_sh": float(np.mean(d <= vsh + 1e-12)),
        "q(sqrt n D) this run": [float(v) for v in np.quantile(sd, QS)],
        "q(sqrt n D) DT0-d f0": [float(v) for v in np.quantile(dt0d, QS)],
        "q(sqrt(0.9) V(B)) m=16000": [float(v) for v in np.quantile(pred, QS)],
        "q(sqrt(0.9) V) with (1-w) instead": [float(v) for v in np.quantile(0.9 * vref, QS)],
        "uniform sqrt n q95 (n=12800 ref)": q95u,
        "rej observed DT0-d": float(np.mean(dt0d >= q95u)),
        "rej predicted sqrt(1-w)": float(np.mean(pred >= q95u)),
        "rej predicted (1-w)": float(np.mean(0.9 * vref >= q95u)),
    }
    out["seconds"] = time.time() - t0
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "path_smoke.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
