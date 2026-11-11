"""Step 14 (DT13) Part N4: Lemma N4 upper bound on the n = 4 SIZE of the uniform-calibrated dip test.

Lemma N4 (DT13-PREREG section 4 v2; proof by chord slopes of the quantile function H = F^{-1}, concave on (0, u],
convex on [u, 1), u = F(m)): with a < b < c < d the order statistics of 4 iid U(0,1) and g1 = b - a, g2 = c - b,
g3 = d - c, for every F in U_c and every r > 1 (FALSE at r = 1),
  P_F(T >= r) <= P(E_u),  E_u = {phi(u) >= r},  phi(v) = A(v)/g1 + (g2 - A(v))/g3,  A(v) = min(max(v - b, 0), g2).
(phi(v) = g2/g3 for v <= b, g2/g1 for v >= c, linear in between: the three cases of the proof.)
At n = 4, {D >= q_4(1-a)} = {T >= r_a} (H&H 1985 spacing reduction; D = T/(4(T+1))).
Reflection (x -> -x) maps u to 1 - u and E_u to E_{1-u}, so u in [1/2, 1] suffices. phi is monotone in v
(constant-linear-constant, slope 1/g1 - 1/g3 on [b, c]), so for a block I = [L, R]: union over v in I of E_v =
E_L u E_R =: E_I.
So size_4(alpha) <= max_j P(E_{I_j}); each P(E_I) is estimated with a one-sided Clopper-Pearson upper limit.
r: the rational r_lo <= r_a from dt13_n4_exact.r_bracket (a smaller r gives a larger event: conservative).
Usage: python dt13_n4_bound.py root N [--quick] -> deposits/dt13/n4_bound_r{root}.json
"""
import json
import os
import sys

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dt13_n4_exact import ALPHAS, r_bracket  # noqa: E402

NBLK = 100           # blocks of width 0.005 on [1/2, 1]
DELTA = 0.001        # simultaneous over 3 alphas x NBLK blocks (one-sided)


def phi(U, v):
    a, b, c, d = U.T
    g1, g2, g3 = b - a, c - b, d - c
    A = np.minimum(np.maximum(v - b, 0.0), g2)
    return A / g1 + (g2 - A) / g3


def counts(U, r, edges):
    """Per block I = [L, R]: #{E_L or E_R} (= #{union over v in I of E_v})."""
    ph = [phi(U, v) >= r for v in edges]
    return [int(np.sum(ph[j] | ph[j + 1])) for j in range(len(edges) - 1)]


def main():
    root, N = int(sys.argv[1]), int(sys.argv[2])
    edges = np.linspace(0.5, 1.0, NBLK + 1)
    rng = np.random.default_rng(np.random.SeedSequence([root, 4]))
    tot = {float(a): np.zeros(NBLK, np.int64) for a in ALPHAS}
    B = 1_000_000
    for s in range(0, N, B):
        m = min(B, N - s)
        U = np.sort(rng.random((m, 4)), axis=1)
        for a in ALPHAS:
            r = float(r_bracket(a)[0])
            tot[float(a)] += np.array(counts(U, r, edges), np.int64)
    conf = 1 - DELTA / (len(ALPHAS) * NBLK)
    res = {"root": root, "N": N, "blocks": edges.tolist(), "delta": DELTA, "alphas": {}}
    for a, k in tot.items():
        up = [1.0 if kk == N else float(stats.beta.ppf(conf, kk + 1, N - kk)) for kk in k]
        j = int(np.argmax(up))
        res["alphas"][str(a)] = {"counts": k.tolist(), "cp_upper": up, "U": up[j], "argmax_block": [edges[j], edges[j + 1]],
                                 "point_max": float(k.max() / N)}
        print(f"alpha={a}: n=4 size <= {up[j]:.5f} (block [{edges[j]:.3f},{edges[j + 1]:.3f}], point {k.max() / N:.5f})")
    od = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "deposits", "dt13")
    try:
        os.makedirs(od, exist_ok=True)
        with open(os.path.join(od, f"n4_bound_r{root}.json"), "w") as fh:
            json.dump(res, fh, indent=1)
    except OSError as err:
        sys.exit(f"cannot write: {err}")


if __name__ == "__main__":
    main()
