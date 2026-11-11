"""DT10 (Step 11, task 2): Monte Carlo certificate of C-edge-infinity on the middle range (DT10-EDGE.md section 6-7).

Target: for each eps-interval I_i = [e_i, e_{i+1}] and level (alpha, t_lo),
    g(1 - eps^2, c*) <= Btilde(eps, c*) <= Btilde(eps, t_lo) <= E[V_i] < alpha   for all eps in I_i,
with Btilde(eps, t) = P(eps G2 > 2t) + E[1{eps G2 <= 2t} Q((2t - eps Lambda(t/eps))/a)], Q = 1 - F_M (Lemma B),
and the pathwise bound (Lemma Lambda0 (ii), (iv); a_i = sqrt(1 - e_i^2) is the larger denominator):
    V_i = 1                                         if e_{i+1} (G_m + 2E) >= 2 t_lo   (possibly infeasible),
    V_i = 1 - F_lo((2 t_lo - e_{i+1} Lam+) / a_i)   otherwise, Lam+ = Lambda_{b_m}(t_lo/e_{i+1} - E) + E.
b_m: piecewise-linear bridge on m equal subintervals; E = ||b - b_m|| drawn exactly (max of m scaled Kolmogorov
suprema, independent of the skeleton). F_lo: Arb lower enclosure of F_M at grid points (floor lookup), from
deposits/dt10/fm_arb_low.npz (fm_cert_check2.py low: grid 0.5..8.5; F_lo = 0 below 0.5).
Usage: python dt10_cert.py m R shard nshards root tag  -> deposits/dt10/cert_<tag>_m<m>_r<root>_s<shard>.json
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dt10_edge_mc import bridge, interp_err  # noqa: E402
from dt10_exact import SCALE, gcm_gap_int, lam_int, to_int  # noqa: E402

ROOT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
K = 36
E_LO, E_HI = 0.12, 0.70710679         # E_HI > 1/sqrt(2) = 0.7071067811865476 (outward rounding)
EDGES = [E_LO * (E_HI / E_LO) ** (i / K) for i in range(K + 1)]
EDGES[-1] = E_HI
LEVELS = ((0.05, 0.7305), (0.01, 0.8585))   # (alpha, rigorous lower end t_lo of c*_alpha; DT8-RESULT)
SAFETY = 1e-9          # on the CDF argument
VAL_PAD = 1e-12        # on every computed G / Lambda value (float error of the exact-predicate kernels <= ~1e-14)
E_ROUND = 2.0 ** -40   # skeleton rounded to the 2^-40 grid (dt10_exact.py)


def load_flo():
    d = np.load(os.path.join(ROOT_DIR, "deposits", "dt10", "fm_arb_low.npz"))
    ys, lo = np.asarray(d["ys"], float), np.asarray(d["lo"], float)
    if not np.all(np.diff(ys) > 0):
        sys.exit("fm_arb_low grid not increasing")
    return ys, np.maximum.accumulate(lo)      # F_M increasing: F(x) >= lo_k for every grid y_k <= x


def f_lo(ys, lo, x):
    k = int(np.searchsorted(ys, x, side="right")) - 1
    return 0.0 if k < 0 else max(0.0, float(lo[k]) - VAL_PAD)    # float(arb.lower()) may round up by an ulp


def path_values(b, E, ys, lo):
    """Conservative V_i for one path. Unresolved cases (integer range exceeded, non-finite) give V = 1."""
    out = np.ones((len(LEVELS), K))
    B, ok = to_int(b)
    if not ok:
        return out, 1
    E = E + E_ROUND
    Gu = gcm_gap_int(B) + VAL_PAD
    for li, (_, t) in enumerate(LEVELS):
        for i in range(K):
            e_lo, e_hi = EDGES[i], EDGES[i + 1]
            if e_hi * (Gu + 2.0 * E) >= 2.0 * t - VAL_PAD:
                continue                                   # V = 1
            y0 = t / e_hi - E
            Y2 = int(np.floor(2.0 * y0 * SCALE)) - 2       # 2y' = Y2/2^40 < 2 y0 (Lambda non-increasing in y)
            if not 0 < Y2 <= 2 ** 44:
                continue
            lam = lam_int(B, Y2) + VAL_PAD + E
            arg = (2.0 * t - e_hi * lam) / np.sqrt(1.0 - e_lo * e_lo) - SAFETY
            v = 1.0 - f_lo(ys, lo, arg)
            out[li, i] = v if np.isfinite(v) else 1.0
    return np.clip(out, 0.0, 1.0), 0


def main():
    m, R, shard, ns = (int(x) for x in sys.argv[1:5])
    root, tag = int(sys.argv[5]), sys.argv[6]
    ys, lo = load_flo()
    reps = R // ns + (1 if shard < R % ns else 0)
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence([root, m, shard])))
    S1 = np.zeros((len(LEVELS), K))
    S2 = np.zeros((len(LEVELS), K))
    t0 = time.time()
    unresolved = 0
    for _ in range(reps):
        b = bridge(rng, m)
        E = interp_err(rng, m)
        v, bad = path_values(b, E, ys, lo)
        unresolved += bad
        S1 += v
        S2 += v * v
    res = {"m": m, "R": R, "shard": shard, "nshards": ns, "root": root, "tag": tag, "reps": reps,
           "unresolved_paths": unresolved,
           "edges": EDGES, "levels": LEVELS, "S1": S1.tolist(), "S2": S2.tolist(), "sec": time.time() - t0}
    out = os.path.join(ROOT_DIR, "deposits", "dt10", f"cert_{tag}_m{m}_r{root}_s{shard}.json")
    try:
        os.makedirs(os.path.dirname(out), exist_ok=True)
        with open(out + ".tmp", "w") as fh:
            json.dump(res, fh)
        os.replace(out + ".tmp", out)
    except OSError as err:
        sys.exit(f"cannot write {out}: {err}")
    print(f"tag={tag} m={m} shard={shard} reps={reps} {res['sec']:.1f}s")


if __name__ == "__main__":
    main()
