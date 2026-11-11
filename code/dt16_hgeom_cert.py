"""Step 16 (DT15) A1': Monte Carlo upper bounds of H_geom(r) = P((1 - v_J) Delta_J < r) by EXACT sampling (DT15-NEPRIME §6).

Law (P6 = Pitman & Uribe Bravo 2012 Thm 1, with Lemma H(a): the bridge has the faces of W on [0, 1], slopes shifted by a
constant, so jumps and order are unchanged): face lengths in size-biased order are uniform stick-breaking
L_1 = U_1, L_i = U_i prod_{j<i}(1 - U_j); given L, slopes s_i = Z_i / sqrt(L_i), Z_i iid N(0, 1); faces are ordered left to
right by increasing slope. J = face 1 (the face containing an independent uniform time). Then
    1 - v_J = sum_{i: s_i > s_1} L_i,   Delta_J = min_{i: s_i > s_1} s_i - s_1   (J last face => 1 - v_J = 0 => A = 0).
Truncation at K faces (rest R_K = prod (1 - U_j)), handled conservatively:
  * ignored faces only ADD length to 1 - v_J  => computed 1 - v_J <= true  (more small A: conservative);
  * an ignored face can make Delta_J smaller only if its slope lies in (s_1, s_1 + Delta_comp). With Delta_comp capped at
    DCAP (count A < r whenever Delta_comp > DCAP), P(that) <= E sum_{i>K} sqrt(L_i / (2 pi)) * DCAP
    = (2/3)^K * 2 / sqrt(2 pi) * DCAP  (E sqrt U = 2/3; self-similarity; E sum_i sqrt L_i = 2), added to every bound.
Per r in RS: n_r = #{A < r}; Clopper-Pearson one-sided upper bound at level DELTA / len(RS) (Bonferroni), plus truncation.
Usage: python dt16_hgeom_cert.py N root tag [DELTA]   -> deposits/dt16/hgeom_<tag>_r<root>.json
"""
import json
import math
import os
import sys
import time

import numpy as np
from scipy.stats import beta

ROOT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
K = 160
DCAP = 1e6
RS = tuple(float(f"{v:.3g}") for v in np.geomspace(1e-5, 0.5, 48))
CHUNK = 50_000


def sample_counts(rng, n):
    u = rng.uniform(size=(n, K))
    rest = np.cumprod(1.0 - u, axis=1)
    L = u.copy()
    L[:, 1:] *= rest[:, :-1]
    s = rng.standard_normal((n, K)) / np.sqrt(np.maximum(L, 1e-300))
    above = s[:, 1:] > s[:, :1]
    one_minus_v = np.sum(np.where(above, L[:, 1:], 0.0), axis=1)
    gap = np.where(above, s[:, 1:] - s[:, :1], np.inf).min(axis=1)
    A = np.where(np.isinf(gap), 0.0, one_minus_v * gap)
    A = np.where(gap > DCAP, -1.0, A)                      # capped: counted as A < r for every r
    return np.array([np.count_nonzero(A < r) for r in RS]), float(rest[:, -1].max())


def main():
    n, root, tag = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    delta = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-4
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([root, 16, 7])))
    cnt, done, rmax, t0 = np.zeros(len(RS), np.int64), 0, 0.0, time.time()
    while done < n:
        m = min(CHUNK, n - done)
        c, r = sample_counts(rng, m)
        cnt += c
        rmax = max(rmax, r)
        done += m
    trunc = (2.0 / 3.0) ** K * 2.0 / math.sqrt(2 * math.pi) * DCAP
    d1 = delta / len(RS)
    ub = [float(beta.ppf(1 - d1, k + 1, n - k)) if k < n else 1.0 for k in cnt.tolist()]
    res = {"N": n, "root": root, "tag": tag, "K": K, "DCAP": DCAP, "delta": delta, "rs": RS, "counts": cnt.tolist(),
           "phat": (cnt / n).tolist(), "cp_upper": ub, "trunc": trunc, "bound": [b + trunc for b in ub],
           "max_rest_RK": rmax, "sec": time.time() - t0}
    out = os.path.join(ROOT_DIR, "deposits", "dt16", f"hgeom_{tag}_r{root}.json")
    try:
        with open(out + ".tmp", "w") as fh:
            json.dump(res, fh)
        os.replace(out + ".tmp", out)
    except OSError as err:
        sys.exit(f"cannot write {out}: {err}")
    for r, p, b in zip(RS, res["phat"], res["bound"], strict=True):
        print(f"r={r:g}: phat {p:.5f}  bound {b:.5f}")
    print(f"N {n} K {K} trunc {trunc:.2e} max R_K {rmax:.2e} {res['sec']:.0f}s -> {out}")


if __name__ == "__main__":
    main()
