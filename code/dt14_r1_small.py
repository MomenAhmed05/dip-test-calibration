"""Step 15 / DT14 WP-D (T15-D): R1 at n = 8..49 via Cor U1' (DT7-UNIF) — certification sampler.

Cor U1' (PROVEN): sup_{F in U_c} P_F(D_n >= d) <= max_j P(Sbar_{n,I_j} >= d), blocks I_j covering [1/2, 1].
Per sample x (sorted iid uniforms) and block I_j we compute R_j = max over leaf pieces [a, b] = I_j cap
[x_(k), x_(k+1)] (x_(0) = 0, x_(n+1) = 1) of piece_ub (certified primal upper bound with proved shape repair,
dt9_block_cert), by branch and bound: a node whose certified bound is <= max(FLOOR, best leaf so far)
is pruned (it cannot raise R_j above that value). If every node is pruned at or below FLOOR we store FLOOR
(meaning "R_j <= FLOOR"), and FLOOR is below every ladder cutoff. LP failure: +inf at inner nodes (split further);
+inf stored at a failing leaf (counts as an exceedance at every cutoff).
Exceedance at cutoff d: R_j >= d - 1e-12 (inclusive; more exceedances can only raise the bound).
Output per chunk: an int16-free compact array of R (CH x J, float64) -> deposits/dt14/r1small/<root>/n<n>_c<k>.npy.
Seeds: chunk k of n uses SeedSequence([root, n, k]); numpy PCG64. Reruns skip existing chunks (resume, same root).
Usage: python dt14_r1_small.py root n chunks [CH]
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dt9_block_cert import piece_ub  # noqa: E402
from sup_pinned import pieces  # noqa: E402

EDGES = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.93, 0.95, 0.97, 0.98, 0.99, 0.995, 1.0]
LADDER_RHO = [round(0.900 + 0.005 * i, 3) for i in range(20)]   # 0.900, 0.905, .., 0.995
C_LO = {0.10: 0.667, 0.05: 0.7313, 0.01: 0.8585}           # lower ends of the certified LF enclosures


def floor_of(n):
    """Strictly below every ladder cutoff at this n."""
    return 0.999 * min(LADDER_RHO) * min(C_LO.values()) / np.sqrt(n)


def block_value(x, P, floor):
    best = floor
    stack = [(0, len(P) - 1)]
    while stack:
        i, j = stack.pop()
        try:
            ub = piece_ub(x, P[i][0], P[j][1])
        except RuntimeError:
            ub = float("inf")
        if ub <= best:
            continue
        if i == j:
            best = ub
            continue
        m = (i + j) // 2
        stack.append((m + 1, j))
        stack.append((i, m))
    return best


def row(x, floor):
    return [block_value(x, pieces(x, lo, hi), floor) for lo, hi in zip(EDGES[:-1], EDGES[1:], strict=True)]


def main():
    root, n, chunks = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    CH = int(sys.argv[4]) if len(sys.argv) > 4 else 5000
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "deposits", "dt14", "r1small", str(root))
    os.makedirs(out, exist_ok=True)
    fl = floor_of(n)
    for k in range(chunks):
        p = os.path.join(out, f"n{n}_c{k}.npy")
        if os.path.exists(p):
            continue
        t0 = time.time()
        rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([root, n, k])))
        R = np.empty((CH, len(EDGES) - 1))
        for i in range(CH):
            R[i] = row(np.sort(rng.random(n)), fl)
        np.save(p + ".tmp.npy", R)
        os.replace(p + ".tmp.npy", p)
        print(f"n={n} chunk {k} {time.time() - t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
