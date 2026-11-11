"""Step 15 / DT14 WP-D at n = 4 (R1): certify that the dip cutoffs d_a = 0.2335 / 0.2420 / 0.2488 (a = 10 / 5 / 1%)
give size <= a at n = 4 uniformly over continuous unimodal laws. Prereg: DT14-R1N4-PREREG.md (frozen at G-029).

Chain (all inputs named):
  P_F(D_4 >= d) = P_F(T >= T_d), T_d = 4d/(1-4d)            [H&H 1985 reduction D_4 = T/(4(T+1)), cited]
               <= sup_{u in [1/2,1]} P(E_u),  E_u = {phi(u) >= T_d}   [Lemma N4 v2, r > 1, PROVEN-written]
               <= max_j P(E_L_j or E_R_j)  over blocks [L_j, R_j] covering [1/2, 1]   [phi monotone in v]
               <= CP_j (one-sided Clopper-Pearson upper limit, simultaneous level 1 - DELTA over 3 x NBLK events).
PASS iff max_j CP_j <= a for each a.  (Order statistics of 4 iid U(0,1); phi as in dt13_n4_bound.py.)
Usage: python dt14_r1_n4.py <root> <N> [--quick]  ->  deposits/dt14/r1n4_r<root>.json
"""
import json
import os
import sys
import time
from fractions import Fraction

import numpy as np
from scipy import stats

CUTS = {0.10: Fraction(2335, 10000), 0.05: Fraction(2420, 10000), 0.01: Fraction(2488, 10000)}
NBLK = 100                 # blocks of width 0.005 on [1/2, 1]
DELTA = 0.001              # simultaneous over 3 x NBLK one-sided limits (Bonferroni)
CHUNK = 2_000_000


def t_of(d):
    return 4 * d / (1 - 4 * d)


REL_TOL = 1e-12            # conservative (inclusive) event classification against float error


def hit(U, v, r):
    """Indicator of E_v = {phi(v) >= r} in the division-free form A g3 + (g2 - A) g1 >= r g1 g3, counted
    INCLUSIVELY within a relative tolerance (more hits => the upper bound stays valid)."""
    a, b, c, d = U.T
    g1, g2, g3 = b - a, c - b, d - c
    A = np.minimum(np.maximum(v - b, 0.0), g2)
    return A * g3 + (g2 - A) * g1 >= r * g1 * g3 * (1.0 - REL_TOL)


def main(root, N, quick):
    t0 = time.time()
    rng = np.random.default_rng(root)
    edges = np.array([0.5 + j / 200 for j in range(NBLK + 1)])
    thr = {a: float(t_of(d)) for a, d in CUTS.items()}
    cnt = {a: np.zeros(NBLK, dtype=np.int64) for a in CUTS}
    done = 0
    while done < N:
        m = min(CHUNK, N - done)
        U = np.sort(rng.random((m, 4)), axis=1)
        for a, T in thr.items():
            h = [hit(U, v, T) for v in edges]
            cnt[a] += np.array([np.sum(h[j] | h[j + 1]) for j in range(NBLK)])
        done += m
    k_events = 3 * NBLK
    lev = 1 - DELTA / k_events
    out = {"root": root, "N": N, "quick": quick, "delta": DELTA, "nblk": NBLK, "rel_tol": REL_TOL,
           "numpy": np.__version__, "bit_generator": type(rng.bit_generator).__name__,
           "edges": "u_j = 1/2 + j/200, j = 0..100", "levels": {}}
    allpass = True
    for a, d in CUTS.items():
        c = cnt[a]
        cp = [float(stats.beta.ppf(lev, int(x) + 1, N - int(x))) if x < N else 1.0 for x in c]
        j = int(np.argmax(cp))
        ok = max(cp) <= a
        allpass &= ok
        out["levels"][str(a)] = {"dip_cut": str(d), "T_cut": str(t_of(d)), "max_phat": float(c.max() / N),
                                 "max_cp_upper": max(cp), "argmax_block": [float(edges[j]), float(edges[j + 1])],
                                 "pass": bool(ok), "counts": [int(x) for x in c], "cp_upper": cp}
    out["result"] = "PASS" if allpass else "FAIL"
    out["seconds"] = round(time.time() - t0, 1)
    here = os.path.dirname(os.path.abspath(__file__))
    tag = "_quick" if quick else ""
    path = os.path.join(here, "..", "deposits", "dt14", f"r1n4_r{root}{tag}.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    q = "--quick" in sys.argv
    main(int(sys.argv[1]), int(sys.argv[2]), q)
