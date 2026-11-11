"""DT9 (Step 11 = next step 1 of the Step 10 report): larger frozen certificate at n = 800 and 1600.

Theorem U1 (DT7-UNIF): sup_F P_F(sqrt(n) D_n > c) <= max_I P(sqrt(n) Sbar_I > c), I over 4 blocks of u.
Here we need only a CONSERVATIVE INDICATOR Z >= 1{sqrt(n) Sbar_I > c} (upper bounds can prove "no",
not "yes"; every sample not certified non-exceeding counts as an exceedance). Decision branch and bound:
  node [P_i .. P_j] of spacing pieces -> certified upper bound ub >= d(G_n, C_[a_i, b_j]) >= every piece inside;
  ub <= c/sqrt(n): no piece inside exceeds (prune); leaf with ub > c/sqrt(n): counted as exceedance (stop).
piece_ub solves the same LP as sup_pinned.piece_dist (same constraints, same solver chain) but returns a
CERTIFIED PRIMAL UPPER BOUND: the band error of the returned knot vector g (both one-sided values
of G_n at every knot), plus the repair cost of its shape violations, plus a rounding allowance, plus SAFETY.
Repair: J_k = slope change of g at knot k. Convex side (knots t_k with k + 1 <= ib, i.e. t_k <= a; LP rows
j + 2 <= ib): delta_k = (-J_k)^+, add delta_k (t_k - s)^+ (convex, zero on [a, 1]). Concave side (knots t_k with
k - 1 >= ia, i.e. t_k >= b): delta'_k = (J_k)^+, add -delta'_k (s - t_k)^+ (concave, zero on [0, b]). The sum is 0
on [a, b], so g + h is in C_[a, b] (continuous; convex on [0, b]; concave on [a, 1]) and |h| <= sum delta + sum
delta'. Rounding allowance 16u(sum |slopes| + K rep + 1), u = 2^-53, covers the floating-point error of the
slope/repair evaluation; SAFETY = 1e-9 covers band and threshold rounding (~1e-17). Non-finite g -> failure.
One-sided rule: a "no exceedance" answer is always backed by a valid upper bound; false "exceedance" answers
only cost power.
LP failure after all fallbacks: at an inner node the bound is +inf (split further); at a leaf the piece is
counted as an exceedance at every cutoff and listed in "failed".
Cutoffs CUTS (sqrt(n) scale) are nested: a higher one is tested only if the lower one is exceeded.
Seeds: chunk k of shard s uses SeedSequence([root, n, s, k]), CH samples per chunk, one JSON per chunk;
a rerun skips chunks whose file exists (resume with the same root; never fresh seeds).
Usage: python dt9_block_cert.py n chunks shard nshards root tag [CH]
"""
import json
import os
import sys
import time

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pinned_dip import V_left  # noqa: E402
from sup_pinned import SAFETY, TRIES, _knots, pieces  # noqa: E402

BLOCKS = [(0.5, 0.9), (0.9, 0.97), (0.97, 0.995), (0.995, 1.0)]
CUTS = (0.7315, 0.86, 0.87)
STATS = {"lps": 0, "fallbacks": 0, "max_repair": 0.0}


U = 2.0 ** -53


def solve(x, a, b):
    """LP of sup_pinned.piece_dist; returns (ts, lo, hi, ia, ib, g). Raises RuntimeError on failure."""
    ts, lo, hi, ia, ib = _knots(x, a, b)
    K = len(ts)
    E = K
    rows, cols, vals, rhs = [], [], [], []

    def add(cs, vs, r):
        k = len(rhs)
        rows.extend([k] * len(cs))
        cols.extend(cs)
        vals.extend(vs)
        rhs.append(r)

    for j in range(K):
        add([j, E], [-1.0, -1.0], -hi[j])
        add([j, E], [1.0, -1.0], lo[j])
    h = np.diff(ts)
    for j in range(K - 2):
        h1, h2 = h[j], h[j + 1]
        if j + 2 <= ib:
            add([j, j + 1, j + 2], [-1.0 / h1, 1.0 / h1 + 1.0 / h2, -1.0 / h2], 0.0)
        if j >= ia:
            add([j, j + 1, j + 2], [1.0 / h1, -1.0 / h1 - 1.0 / h2, 1.0 / h2], 0.0)
    A = coo_matrix((vals, (rows, cols)), shape=(len(rhs), E + 1)).tocsr()
    c = np.zeros(E + 1)
    c[E] = 1.0
    bnds = [(None, None)] * E + [(0, None)]
    res = None
    STATS["lps"] += 1
    for meth, opts in TRIES:
        res = linprog(c, A_ub=A, b_ub=np.array(rhs), bounds=bnds, method=meth, options=opts)
        if res.status == 0:
            break
        STATS["fallbacks"] += 1
    if res is None or res.status != 0:
        raise RuntimeError("LP failed after fallbacks")
    g = res.x[:E]
    if not np.all(np.isfinite(g)):
        raise RuntimeError("non-finite LP solution")
    return ts, lo, hi, ia, ib, g


def bound_of(ts, lo, hi, ia, ib, g):
    """Certified upper bound for a returned knot vector g (see module docstring)."""
    K = len(ts)
    band = float(np.max(np.maximum(hi - g, g - lo)))
    sl = np.diff(g) / np.diff(ts)
    ds = np.diff(sl)                      # ds[j] = slope change at knot j + 1
    jj = np.arange(K - 2)
    rep = float(np.sum(np.maximum(-ds[jj + 2 <= ib], 0.0))) + float(np.sum(np.maximum(ds[jj >= ia], 0.0)))
    rnd = 16.0 * U * (float(np.sum(np.abs(sl))) + K * rep + 1.0)
    ub = band + rep + rnd + SAFETY
    if not np.isfinite(ub):
        raise RuntimeError("non-finite bound")
    STATS["max_repair"] = max(STATS["max_repair"], rep)
    return ub


def piece_ub(x, a, b):
    """Certified primal upper bound on d(G_n, C_[a,b]) (raw scale). Raises RuntimeError on failure."""
    return bound_of(*solve(x, a, b))


def exceeds(x, P, t):
    """Conservative decision: returns False only if a certified upper bound shows every piece in P has
    d(G_n, C_piece) <= t; otherwise True (exceedance or unresolved). Second value: an LP failed."""
    failed = False
    stack = [(0, len(P) - 1)]
    while stack:
        i, j = stack.pop()
        try:
            ub = piece_ub(x, P[i][0], P[j][1])
        except RuntimeError:
            if i == j:
                return True, True
            ub = float("inf")
            failed = True
        if ub <= t:
            continue
        if i == j:
            return True, failed
        m = (i + j) // 2
        stack.append((m + 1, j))
        stack.append((i, m))
    return False, failed


def sample_row(x):
    """Per block: list of exceedance indicators for CUTS (nested), and a failure flag."""
    rt = np.sqrt(len(x))
    row, fl = [], []
    for lo_u, hi_u in BLOCKS:
        P = pieces(x, lo_u, hi_u)
        ex, anyf = [], False
        for c in CUTS:
            if ex and not ex[-1]:
                ex.append(False)
                continue
            e, f = exceeds(x, P, c / rt)
            anyf = anyf or f
            ex.append(bool(e))
        row.append(ex)
        fl.append(anyf)
    return row, fl


def main():
    n, chunks, shard, ns = (int(a) for a in sys.argv[1:5])
    root = int(sys.argv[5])
    tag = sys.argv[6]
    CH = int(sys.argv[7]) if len(sys.argv) > 7 else 500
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "deposits", "dt9")
    try:
        os.makedirs(out, exist_ok=True)
    except OSError as err:
        sys.exit(f"cannot create {out}: {err}")
    rt = np.sqrt(n)
    for k in range(chunks):
        path = os.path.join(out, f"{tag}_n{n}_r{root}_s{shard}_c{k:04d}.json")
        if os.path.exists(path):
            continue
        rng = np.random.default_rng(np.random.SeedSequence([root, n, shard, k]))
        t0 = time.time()
        for key in ("lps", "fallbacks"):
            STATS[key] = 0
        STATS["max_repair"] = 0.0
        R, F, V = [], [], []
        for _ in range(CH):
            x = np.sort(rng.uniform(size=n))
            r, f = sample_row(x)
            R.append(r)
            F.append(f)
            V.append(rt * V_left(x, 1.0))
        res = {"n": n, "root": root, "shard": shard, "nshards": ns, "chunk": k, "CH": CH, "blocks": BLOCKS,
               "cuts": CUTS, "ex": R, "failed": F, "V": V, "stats": dict(STATS), "seconds": time.time() - t0}
        tmp = path + ".tmp"
        try:
            with open(tmp, "w") as fh:
                json.dump(res, fh)
            os.replace(tmp, path)
        except OSError as err:
            sys.exit(f"cannot write {path}: {err}")
        print(f"n={n} shard={shard} chunk={k} {time.time() - t0:.1f}s lps={STATS['lps']} "
              f"fb={STATS['fallbacks']} rep={STATS['max_repair']:.1e}", flush=True)


if __name__ == "__main__":
    main()
