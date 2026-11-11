"""Step 15 / DT14 WP-D (T15-D) analysis of dt14_r1_small.py output (frozen rule, DT14-R1SMALL-PREREG).

Per (n, alpha): ladder d_l = rho_l * C_LO[alpha] / sqrt(n), rho_l increasing. For each l and block j:
K_lj = #{R_j >= d_l - 1e-12}, one-sided Clopper-Pearson upper U_lj = BetaInv(1 - eta; K + 1, N - K) (1 if K = N),
eta = DELTA / (#n * 3). PASS at l iff max_j U_lj <= alpha. Reported cutoff = first passing l (fixed-sequence
argument, PREREG §3). "Can reject" iff d < dip of the explicit two-cluster example: 1/4 (n even), (n-1)/(4n) (n odd).
Usage: python dt14_r1_small_analysis.py root N out.json
"""
import glob
import json
import os
import sys

import numpy as np
from scipy.stats import beta

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dt14_r1_small import C_LO, LADDER_RHO  # noqa: E402

NS = list(range(8, 50))
DELTA = 0.001


def cp_upper(K, N, eta):
    return 1.0 if K >= N else float(beta.ppf(1.0 - eta, K + 1, N - K))


def two_cluster_dip(n):
    return 0.25 if n % 2 == 0 else (n - 1) / (4.0 * n)


def main():
    root, N, outp = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "deposits", "dt14", "r1small", str(root))
    eta = DELTA / (len(NS) * 3)
    res = {"root": root, "N_per_n": N, "eta": eta, "delta": DELTA, "cells": {}}
    for n in NS:
        files = sorted(glob.glob(os.path.join(base, f"n{n}_c*.npy")), key=lambda s: int(s.rsplit("_c", 1)[1][:-4]))
        R = np.concatenate([np.load(f) for f in files]) if files else np.empty((0, 15))
        if len(R) < N:
            res["cells"][str(n)] = {"status": f"INCOMPLETE {len(R)}/{N}"}
            continue
        R = R[:N]
        cell = {"n_inf": int(np.sum(~np.isfinite(R)))}
        for a, c in C_LO.items():
            ladder = []
            first = None
            for rho in LADDER_RHO:
                d = rho * c / np.sqrt(n)
                K = (R >= d - 1e-12).sum(axis=0)
                U = [cp_upper(int(k), N, eta) for k in K]
                ok = max(U) <= a
                ladder.append({"rho": rho, "d": d, "maxK": int(K.max()), "argmax_block": int(np.argmax(U)),
                               "maxU": max(U), "pass": ok})
                if ok and first is None:
                    first = {"rho": rho, "d": d, "maxU": max(U), "phat": float(K.max() / N),
                             "lf_lo": c / np.sqrt(n), "can_reject": bool(d < two_cluster_dip(n))}
            cell[str(a)] = {"certified": first, "ladder": ladder}
        res["cells"][str(n)] = cell
    try:
        with open(outp, "w") as fh:
            json.dump(res, fh, indent=1)
    except OSError as err:
        sys.exit(f"cannot write {outp}: {err}")
    for n, cell in res["cells"].items():
        if "status" in cell:
            print(n, cell["status"])
            continue
        s = []
        for a in ("0.1", "0.05", "0.01"):
            f = cell[a]["certified"]
            s.append("none" if f is None else f"{f['d']:.4f}(rho{f['rho']:.3f} U{f['maxU']:.4f}{'' if f['can_reject'] else ' NOREJ'})")
        print(f"n={n} inf={cell['n_inf']} " + " ".join(s))


if __name__ == "__main__":
    main()
