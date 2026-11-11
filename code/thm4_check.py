"""Numerical cross-check of DT3-THEOREM §4c Theorem 4 (directional upper bound for the dip).

For a fixed unimodal CDF F on [0, 1] and a fixed smooth h with h(0) = h(1) = 0, compute
  ratio(r) = D(F + r h) / r      (grid LP: best convex-then-concave piecewise-linear fit, free turn node)
and compare it with
  T_F(h) = max( max_{affine J left} V_J(h), D_P(h), max_{affine J right} C_J(h) ).
Theorem 4 says limsup ratio(r) <= T_F(h). For step-like stretches Theorem 3 also gives the matching
lower bound, so we expect ratio(r) -> T_F(h) in every case here.

Deterministic (no RNG). Grid discretisation: values are sampled at the nodes and g is piecewise linear, so
small differences (order of curvature x grid^2 / r) are expected. DEV check, not evidence for a claim tier.
Run: python code/thm4_check.py [--quick]
"""
import sys
import time

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, vstack

QUICK = "--quick" in sys.argv
N = 240 if QUICK else 600
RS = (0.1, 0.03, 0.01) if QUICK else (0.1, 0.03, 0.01, 0.003)


def cdf_from_density(dens, x):
    """CDF on grid x from a vectorised density (trapezoid rule on a fine grid, normalised to F(1) = 1)."""
    fine = np.linspace(0.0, 1.0, 200_001)
    fv = dens(fine)
    cum = np.concatenate([[0.0], np.cumsum((fv[1:] + fv[:-1]) / 2 * np.diff(fine))])
    cum /= cum[-1]
    return np.interp(x, fine, cum)


def fit_lp(x, y, k_range):
    """min over turn node k in k_range of  min_g max_i |y_i - g_i|, g convex on nodes <= k, concave on >= k."""
    n = len(x)
    best = (np.inf, None)
    dx = np.diff(x)
    for k in k_range:
        rows, cols, vals = [], [], []
        r = 0
        # second-difference constraints: sign * (slope_{i} - slope_{i-1}) <= 0 written as A g <= 0
        for i in range(1, n - 1):
            if i == k:
                continue
            sgn = -1.0 if i < k else 1.0  # convex: slope_i - slope_{i-1} >= 0  ->  -(...) <= 0
            a, b = 1.0 / dx[i - 1], 1.0 / dx[i]
            for c, v in ((i - 1, a), (i, -a - b), (i + 1, b)):
                rows.append(r)
                cols.append(c)
                vals.append(sgn * v)
            r += 1
        A_shape = coo_matrix((vals, (rows, cols)), shape=(r, n + 1))
        # |y - g| <= t :  g - t <= y ;  -g - t <= -y
        eye = coo_matrix((np.ones(n), (np.arange(n), np.arange(n))), shape=(n, n + 1))
        tcol = coo_matrix((np.ones(n), (np.arange(n), np.full(n, n))), shape=(n, n + 1))
        A = vstack([A_shape, eye - tcol, -eye - tcol]).tocsr()
        bvec = np.concatenate([np.zeros(r), y, -y])
        cvec = np.zeros(n + 1)
        cvec[n] = 1.0
        res = linprog(cvec, A_ub=A, b_ub=bvec, bounds=[(None, None)] * n + [(0, None)], method="highs")
        if not res.success:
            raise RuntimeError(f"LP failed at k={k}: {res.message}")
        if res.fun < best[0]:
            best = (res.fun, k)
    return best


def half_gap(x, y, convex=True):
    """V (convex=True) or C (convex=False) on the given nodes: half the max gap to the lower/upper hull."""
    s = 1.0 if convex else -1.0
    yy = s * y
    hull = [0]
    for i in range(1, len(x)):
        while len(hull) >= 2:
            i0, i1 = hull[-2], hull[-1]
            if (yy[i1] - yy[i0]) * (x[i] - x[i0]) >= (yy[i] - yy[i0]) * (x[i1] - x[i0]):
                hull.pop()
            else:
                break
        hull.append(i)
    g = np.interp(x, x[hull], yy[hull])
    return 0.5 * float(np.max(yy - g))


def T_F(x, h, left, plateau, right):
    terms = [half_gap(x[(x >= a) & (x <= b)], h[(x >= a) & (x <= b)], True) for a, b in left]
    terms += [half_gap(x[(x >= a) & (x <= b)], h[(x >= a) & (x <= b)], False) for a, b in right]
    if plateau is not None:
        m = (x >= plateau[0]) & (x <= plateau[1])
        xs, hs = x[m], h[m]
        terms.append(fit_lp(xs, hs, range(len(xs)))[0])
    return max(terms) if terms else 0.0


def h_test(x):
    return 0.12 * np.sin(2 * np.pi * x) + 0.08 * np.sin(6 * np.pi * x) + 0.05 * np.sin(11 * np.pi * x)


def h_family(x):
    """h_test plus 6 deterministic sine series (coefficients cos(7j + s)/j, j = 1..12) and a plateau bump."""
    hs = [("h_test", h_test(x))]
    for s in range(6):
        hs.append((f"sine{s}", sum(np.cos(7 * j + s) / j * 0.06 * np.sin(j * np.pi * x) for j in range(1, 13))))
    hs.append(("bump0.67", 0.1 * np.exp(-((x - 0.67) / 0.03) ** 2) * np.sin(40 * x) * x * (1 - x) * 4))
    return hs


# Laws: (name, case, density, affine stretches left, modal plateau or None, affine stretches right, mode, window)
def lin(a, b, ya, yb):
    return lambda t: ya + (yb - ya) * (t - a) / (b - a)


LAWS = [
    ("A: flat shoulder + ramp + plateau + smooth right", "A",
     lambda t: np.where(t < 0.4, 0.5, np.where(t < 0.6, lin(0.4, 0.6, 0.5, 1.8)(t),
                np.where(t < 0.75, 1.8, 1.8 * (1 - (t - 0.75) / 0.25) ** 2))),
     [(0.0, 0.4)], (0.6, 0.75), [], None),
    ("B: smooth rise, jump up at mode 0.55, flat right block", "B",
     lambda t: np.where(t < 0.55, 0.3 + 1.2 * t, np.where(t < 0.8, lin(0.55, 0.8, 1.6, 1.0)(t), 0.6)),
     [], None, [(0.8, 1.0)], 0.55),
    ("C': Beta(2,2) shape, no flat stretch (T = 0)", "C'",
     lambda t: t * (1 - t), [], None, [], 0.5),
    ("C': flat left block + smooth point mode at 0.6", "C'",
     lambda t: np.where(t < 0.35, 0.4, np.where(t < 0.6, 0.4 + 2.0 * (t - 0.35), 0.9 - 0.9 * (t - 0.6))),
     [(0.0, 0.35)], None, [], 0.6),
]


def main():
    t0 = time.time()
    x = np.linspace(0.0, 1.0, N + 1)
    hs = h_family(x)[:2] if QUICK else h_family(x)
    worst = []
    for name, case, dens, left, plateau, right, mode in LAWS:
        F = cdf_from_density(dens, x)
        if plateau is not None:
            lo, hi = plateau
        else:
            lo = hi = mode
        ks = [i for i in range(N + 1) if lo - 0.05 <= x[i] <= hi + 0.05]
        print(f"{name} [case {case}]")
        for hname, h in hs:
            h = h - h[0] - (h[-1] - h[0]) * x  # enforce h(0) = h(1) = 0
            T = T_F(x, h, left, plateau, right)
            ratios = [fit_lp(x, F + r * h, ks)[0] / r for r in RS]
            worst.append((name[:3], hname, ratios[-1] - T))
            print(f"   {hname:<9} T={T:.5f}  D/r: " + "  ".join(f"{q:.5f}" for q in ratios)
                  + f"   last-T={ratios[-1] - T:+.5f}")
    print("max over (law, h) of D/r - T at the smallest r:", f"{max(w[2] for w in worst):+.5f}",
          "| min:", f"{min(w[2] for w in worst):+.5f}", f"| r = {RS}")
    print(f"grid N = {N}, {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
