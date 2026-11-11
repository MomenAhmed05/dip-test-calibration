"""Summarise dt16_pair_cert.py shards: per level and cell, mean Y, SE, empirical-Bernstein upper bound.

Bound (Maurer-Pontil, Y in [-1, 1], range 2): mean + sqrt(2 V log(2/d) / N) + 14 log(2/d) / (3 (N - 1)),
d = DELTA / (number of cells over all levels) (Bonferroni). A cell is certified iff the bound is < 0.
Reach = smallest eps_lo such that every cell from it up to 1/sqrt 2 is certified.
Sign test (adopted): P_C = 2^-K sum_{j <= n_plus} C(K, j), K = n_plus + n_minus, exact integers; K = 0 => 1.
Intersection-union (all-cells-pass) at level DELTA: certified iff max_C P_C < DELTA (no Bonferroni).
Usage: python dt16_pair_cert_summary.py TAG ROOT [DELTA]
"""
import glob
import json
import math
import os
import sys
from fractions import Fraction

import numpy as np

ROOT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def sign_p(n_plus, n_minus):
    """Exact lower-tail sign-test p-value as a float (decision uses exact integers in sign_reject)."""
    k = n_plus + n_minus
    if k == 0:
        return 1.0
    return float(Fraction(sum(math.comb(k, j) for j in range(n_plus + 1)), 2 ** k))


def sign_reject(n_plus, n_minus, delta_num, delta_den):
    """Exact: P_C < delta_num / delta_den with P_C = A / 2^K."""
    k = n_plus + n_minus
    if k == 0:
        return False
    a = sum(math.comb(k, j) for j in range(n_plus + 1))
    return a * delta_den < delta_num * 2 ** k


def main():
    tag, root = sys.argv[1], sys.argv[2]
    delta = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-3
    fs = sorted(glob.glob(os.path.join(ROOT_DIR, "deposits", "dt16", "pair", f"{tag}_r{root}_s*.json")))
    if not fs:
        sys.exit("no shards")
    ds = []
    for f in fs:
        try:
            with open(f) as fh:
                ds.append(json.load(fh))
        except (OSError, ValueError) as err:
            sys.exit(f"bad shard {f}: {err}")
    levels, edges = ds[0]["levels"], ds[0]["edges"]
    n = sum(d["reps"] for d in ds)
    keys = ("m0", "mmax", "R", "nshards", "root", "tag", "levels", "edges")
    complete = (len(ds) == ds[0]["nshards"] and n == ds[0]["R"]
                and sorted(d["shard"] for d in ds) == list(range(ds[0]["nshards"]))
                and all(all(d[k] == ds[0][k] for k in keys) for d in ds))
    print(f"COMPLETE {complete} (shards {len(ds)}/{ds[0]['nshards']}, N {n}/{ds[0]['R']})")
    ncell = sum(len(e) - 1 for e in edges)
    lg = math.log(2.0 / (delta / ncell))
    depths, sec = {}, 0.0
    for d in ds:
        sec += d["sec"]
        for k, v in d["depths"].items():
            depths[int(k)] = depths.get(int(k), 0) + v
    print(f"{tag} root {root}: {len(ds)} shards, N = {n}, cells {ncell}, delta {delta}, CPU {sec:.0f}s "
          f"({sec / n:.2f}s/path), unresolved {sum(d['unresolved_paths'] for d in ds)}, depths {dict(sorted(depths.items()))}")
    for li, (lv, e) in enumerate(zip(levels, edges, strict=True)):
        s1 = np.sum([d["S1"][li] for d in ds], axis=0)
        s2 = np.sum([d["S2"][li] for d in ds], axis=0)
        npl = np.sum([d["n_plus"][li] for d in ds], axis=0)
        nmi = np.sum([d["n_minus"][li] for d in ds], axis=0)
        mean = s1 / n
        var = np.maximum(0.0, (s2 - n * mean ** 2) / (n - 1))
        ub = mean + np.sqrt(2 * var * lg / n) + 14 * lg / (3 * (n - 1))
        ok = ub < 0
        reach = math.nan
        for i in range(len(ok) - 1, -1, -1):
            if not ok[i]:
                break
            reach = e[i]
        num, den = (int(round(delta * 10 ** 9)), 10 ** 9)
        sg = np.array([sign_reject(int(a), int(b), num, den) for a, b in zip(npl, nmi, strict=True)])
        pc = np.array([sign_p(int(a), int(b)) for a, b in zip(npl, nmi, strict=True)])
        sreach = math.nan
        for i in range(len(sg) - 1, -1, -1):
            if not sg[i]:
                break
            sreach = e[i]
        print(f"DECISION level alpha={lv[0]}: {'CERTIFIED' if complete and bool(sg.all()) else 'NOT CERTIFIED'} "
              f"(sign-IUT, every cell P_C < {delta}; cells {len(sg)}; eps [{e[0]}, {e[-1]}])")
        print(f"level alpha={lv[0]} t'={lv[1]} eps_lo={lv[2]}: Bernstein(Bonf) {ok.sum()}/{len(ok)} reach {reach}; "
              f"sign-IUT {sg.sum()}/{len(sg)} reach {sreach} max P_C {pc.max():.3g} (cell {int(pc.argmax())}, "
              f"eps {e[int(pc.argmax())]:.4f}) min K {int((npl + nmi).min())} max n+ {int(npl.max())}")
        for i in sorted(set(np.linspace(0, len(ok) - 1, 12).astype(int).tolist())):
            print(f"  eps[{e[i]:.4f},{e[i + 1]:.4f}] mean {mean[i]:+.5f} SE {math.sqrt(var[i] / n):.5f} "
                  f"P+ {npl[i] / n:.5f} P- {nmi[i] / n:.5f} UB {ub[i]:+.5f} P_C {pc[i]:.3g}")


if __name__ == "__main__":
    main()
