"""O5-ML secondary (GOAL v1.2): how often does the gap statistic (Tibshirani, Walther & Hastie 2001) choose k > 1 on
one-mode pile-up laws vs log-concave controls? 1-D data, ClustPy GapStatistic (k-means, uniform reference box,
'original' 1-SE stopping rule), B = 100 reference sets, k in 1..6. Not a dip method: outside the calibration comparison.
Run (venv): python code/o5_gap.py [--quick] [--shard i/k]
Roots: CONFIRM 2026110665, DEV quick 2026110666.
"""
# pyright: reportMissingImports=false
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from o5_ml import sample_1d  # noqa: E402

QUICK = "--quick" in sys.argv
ROOT = 2026110666 if QUICK else 2026110665
OUT = HERE.parent / "deposits" / ("o5gap_quick" if QUICK else "o5gap")
R = 40 if QUICK else 200  # 1000 -> 200 at freeze (~9 s per fit with B = 100; secondary item)
NS = (100, 500)
LAWS = [("N1", "pileup", ("W", 0.1)), ("N2_s0.002", "pileup", ("clipjit", 0.002)),
        ("R0_s0.002", "pileup", ("relujit", 0.002)), ("A_N2", "pileup-atom", ("clipjit", 0.0)),
        ("U", "logconcave", ("uniform", None)), ("G", "logconcave", ("normal", None)),
        ("logistic", "logconcave", ("logistic", None))]


def draw(kind, par, n, g):
    if kind == "logistic":
        return g.logistic(size=n)
    return sample_1d(kind, par, n, g)


def run_cell(li, n):
    from clustpy.centroid import GapStatistic
    name, role, (kind, par) = LAWS[li]
    k = []
    for r in range(R):
        g = np.random.default_rng(np.random.SeedSequence([ROOT, li, n, r]))
        x = draw(kind, par, n, g).reshape(-1, 1)
        seed = int(np.random.SeedSequence([ROOT, 99, li, n, r]).generate_state(1)[0] % (2 ** 31 - 1))
        k.append(GapStatistic(max_n_clusters=6, n_boots=100, random_state=seed).fit(x).n_clusters_)
    return np.asarray(k)


def main():
    (OUT / "cells").mkdir(parents=True, exist_ok=True)
    shard = sys.argv[sys.argv.index("--shard") + 1] if "--shard" in sys.argv else None
    cells = [(li, n) for li in range(len(LAWS)) for n in NS]
    if shard is not None:
        i, k = (int(v) for v in shard.split("/"))
        t0 = time.time()
        for j, (li, n) in enumerate(cells):
            p = OUT / "cells" / f"L{li}_n{n}.npy"
            if j % k == i and not p.exists():
                np.save(p, run_cell(li, n))
                print(f"[{time.time() - t0:6.0f}s] {LAWS[li][0]} n={n}", flush=True)
        return
    rows = []
    for li, n in cells:
        p = OUT / "cells" / f"L{li}_n{n}.npy"
        if not p.exists():
            print(f"missing {p.name}; run shards first")
            return
        k = np.load(p)
        rows.append({"law": LAWS[li][0], "role": LAWS[li][1], "n": n, "R": int(k.size),
                     "P_k_gt_1": float(np.mean(k > 1)), "k_hist": np.bincount(k, minlength=7).tolist()})
        print(f"{LAWS[li][0]:10s} {LAWS[li][1]:11s} n={n:4d} P(k>1) = {100 * np.mean(k > 1):5.1f}%  hist {np.bincount(k, minlength=7)[1:].tolist()}")
    (OUT / "summary.json").write_text(json.dumps({"root": ROOT, "rows": rows}, indent=1))


if __name__ == "__main__":
    main()
