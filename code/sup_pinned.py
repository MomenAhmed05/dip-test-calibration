"""Block envelope of the pinned-mode dip (Step 10 / DT8; DT7-UNIF.md).

For a sorted uniform sample x on [0, 1] with empirical CDF G_n, and a closed spacing piece
[a, b] (a <= b, [a, b] inside one closed spacing [x_(k), x_(k+1)], x_(0) = 0, x_(n+1) = 1):
  C_[a,b] = {g continuous on [0,1]: convex on [0, b], concave on [a, 1]}  (subset of K_u for all u in [a, b])
  d(G_n, C_[a,b]) = LP over the values of g at the knots {0, x_i, a, b, 1} (piecewise-linear reduction;
  the error is checked against both one-sided values of G_n at every sample knot).
Block statistic for a block I = [lo, hi] of pinned modes:
  Sbar_I = max over spacing pieces [a, b] = I cap [x_(k), x_(k+1)] (nonempty) of d(G_n, C_[a,b]).
Theorem (DT7-UNIF Thm U1): sup_F P_F(sqrt(n) D_n > t) <= max_I P(sqrt(n) Sbar_I > t) for every n.
Numerics: HiGHS LP (minimisation; a feasible-side safety margin is added by the caller, see SAFETY).
"""
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

SAFETY = 1e-9   # added to every LP value (raw scale) before use; LP tolerances are 1e-10
MAXVIOL = [0.0]  # running max of (slope-violation x smallest spacing): must stay << SAFETY
FALLBACKS = [0]  # count of solver retries
TRIES = [("highs", {"primal_feasibility_tolerance": 1e-10, "dual_feasibility_tolerance": 1e-10}),
         ("highs", {}), ("highs-ipm", {}), ("highs-ds", {})]


def _knots(x, a, b):
    """Sorted unique knots with (lo, hi) one-sided values of G_n; indices of a and b."""
    n = len(x)
    ts = np.concatenate([[0.0], x, [1.0]])
    hi = np.concatenate([[0.0], np.arange(1, n + 1) / n, [1.0]])
    lo = np.concatenate([[0.0], np.arange(0, n) / n, [1.0]])
    extra = []
    for p in (a, b):
        if not np.any(ts == p):
            k = int(np.searchsorted(x, p, side="right"))
            extra.append((p, k / n))
    if extra:
        et = np.array([e[0] for e in extra])
        ev = np.array([e[1] for e in extra])
        ts = np.concatenate([ts, et])
        hi = np.concatenate([hi, ev])
        lo = np.concatenate([lo, ev])
        o = np.argsort(ts, kind="stable")
        ts, hi, lo = ts[o], hi[o], lo[o]
        keep = np.concatenate([[True], np.diff(ts) > 0])
        ts, hi, lo = ts[keep], hi[keep], lo[keep]
    ia = int(np.searchsorted(ts, a))
    ib = int(np.searchsorted(ts, b))
    return ts, lo, hi, ia, ib


def piece_dist(x, a, b):
    """d(G_n, C_[a,b]) by LP. x sorted in (0, 1)."""
    ts, lo, hi, ia, ib = _knots(x, a, b)
    K = len(ts)
    E = K
    rows, cols, vals, rhs = [], [], [], []

    def add(cs, vs, b):
        r = len(rhs)
        rows.extend([r] * len(cs))
        cols.extend(cs)
        vals.extend(vs)
        rhs.append(b)

    for j in range(K):                 # hi_j - e <= g_j <= lo_j + e
        add([j, E], [-1.0, -1.0], -hi[j])
        add([j, E], [1.0, -1.0], lo[j])
    h = np.diff(ts)
    for j in range(K - 2):
        h1, h2 = h[j], h[j + 1]
        if j + 2 <= ib:                # convex on [0, b]: s_j <= s_{j+1}
            add([j, j + 1, j + 2], [-1.0 / h1, 1.0 / h1 + 1.0 / h2, -1.0 / h2], 0.0)
        if j >= ia:                    # concave on [a, 1]: s_j >= s_{j+1}
            add([j, j + 1, j + 2], [1.0 / h1, -1.0 / h1 - 1.0 / h2, 1.0 / h2], 0.0)
    A = coo_matrix((vals, (rows, cols)), shape=(len(rhs), E + 1)).tocsr()
    c = np.zeros(E + 1)
    c[E] = 1.0
    bnds = [(None, None)] * E + [(0, None)]
    res = None
    for meth, opts in TRIES:     # tight tolerances first; fall back on solver failure (logged)
        res = linprog(c, A_ub=A, b_ub=np.array(rhs), bounds=bnds, method=meth, options=opts)
        if res.status == 0:
            break
        FALLBACKS[0] += 1
    if res is None or res.status != 0:
        raise RuntimeError(f"LP failed after fallbacks: {res.message if res else ''}")
    g, e = res.x[:E], float(res.x[E])
    # a-posteriori checks: realised band error of the returned g, and its worst shape violation
    band = float(np.max(np.maximum(hi - g, g - lo)))
    s = np.diff(g) / h
    ds = np.diff(s)
    jj = np.arange(K - 2)
    viol = max(float(np.max(-ds[jj + 2 <= ib], initial=0.0)), float(np.max(ds[jj >= ia], initial=0.0)))
    MAXVIOL[0] = max(MAXVIOL[0], viol * float(np.min(h)))   # value-scale size of the worst kink error
    return max(e, band) + SAFETY


def pieces(x, lo_u, hi_u):
    """Closed spacing pieces [a, b] = [lo_u, hi_u] cap [x_(k), x_(k+1)], nonempty."""
    pts = np.concatenate([[0.0], x, [1.0]])
    out = []
    for k in range(len(pts) - 1):
        a, b = max(pts[k], lo_u), min(pts[k + 1], hi_u)
        if a <= b:
            out.append((a, b))
    return out


def block_stat(x, lo_u, hi_u):
    return max(piece_dist(x, a, b) for a, b in pieces(x, lo_u, hi_u))


def block_stat_bb(x, lo_u, hi_u, floor):
    """Branch and bound. Returns (value, exact): exact max over pieces if it is >= floor; otherwise
    an upper bound < floor (exact=False). Uses C_[a_i, b_j] subset of C_[a, b] for [a, b] inside
    [a_i, b_j], so d(G_n, C_[a_i, b_j]) bounds every piece inside from above."""
    P = pieces(x, lo_u, hi_u)
    best = -1.0
    stack = [(0, len(P) - 1, None)]
    ub_all = 0.0
    while stack:
        i, j, ub = stack.pop()
        if ub is None:
            ub = piece_dist(x, P[i][0], P[j][1])
        if ub <= max(best, floor) and ub < floor:
            ub_all = max(ub_all, ub)
            continue
        if ub <= best:
            continue
        if i == j:
            best = max(best, ub)
            continue
        m = (i + j) // 2
        stack.append((i, m, None))
        stack.append((m + 1, j, None))
    if best >= floor:
        return best, True
    return max(best, ub_all), False


if __name__ == "__main__":
    import os
    import sys
    import time
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from pinned_dip import V_left, Y_joined
    rng = np.random.default_rng(2026110702)   # DEV smoke root (CONVENTIONS: DT8)
    worst_lb, worst_edge = 0.0, 0.0
    for n in (30, 60):
        for _ in range(10):
            x = np.sort(rng.uniform(size=n))
            for (a, b) in pieces(x, 0.5, 1.0):
                d = piece_dist(x, a, b)
                for u in np.linspace(a, b, 3):
                    if 0.0 < u < 1.0:
                        worst_lb = max(worst_lb, Y_joined(x, u) - d)
            # last piece [x_(n), 1] at the edge: C = convex on [0,1], linear on [x_(n),1]; G_n = 1 there,
            # so it equals the plain convex distance V
            worst_edge = max(worst_edge, abs(piece_dist(x, x[-1], 1.0) - V_left(x, 1.0)))
    print(f"max (Y(u) - d(C_[a,b])) over u in piece (must be <= 0): {worst_lb:.3e}")
    print(f"max |d(C_[x_n,1]) - V| (should be ~0): {worst_edge:.3e}")
    for n in (200, 400, 800):
        x = np.sort(rng.uniform(size=n))
        print(f"n={n}: pieces in [.5,1]: {len(pieces(x, .5, 1))}")
    # branch-and-bound agrees with brute force (exact part) and its timing
    bad = 0
    for n in (40, 80):
        for _ in range(6):
            x = np.sort(rng.uniform(size=n))
            for (lo_u, hi_u) in ((0.5, 0.8), (0.8, 1.0)):
                bf = block_stat(x, lo_u, hi_u)
                for fl in (0.0, bf - 1e-6, bf + 1e-6):
                    v, ex = block_stat_bb(x, lo_u, hi_u, fl)
                    if (ex and abs(v - bf) > 1e-10) or (not ex and not (bf - 1e-10 <= v < fl)):
                        bad += 1
    print("bb vs brute-force mismatches:", bad)
    for n in (200, 800):
        t0 = time.time()
        for _ in range(3):
            x = np.sort(rng.uniform(size=n))
            for blk in ((0.5, 0.9), (0.9, 0.99), (0.99, 1.0)):
                block_stat_bb(x, *blk, 0.70 / np.sqrt(n))
        print(f"n={n}: bb 3 blocks, mean {(time.time() - t0) / 3:.2f} s/sample")
    print(f"max value-scale shape violation seen: {MAXVIOL[0]:.2e} (SAFETY {SAFETY:.0e})")
