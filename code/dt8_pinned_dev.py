"""DT8 DEV diagnostic (Step 10): tails of the pinned-mode envelope.

For uniform samples of size n, at pinned modes u (symmetry: u and 1-u have the same law), estimate
P(sqrt(n) Y(u) > t) (joined, = the domination-lemma bound) and P(sqrt(n) Z(u) > t) (decoupled), at
t = 0.73139 (LF 5% cutoff) and 0.86 (LF 1%), plus the edge u = 1. Common samples across u.
DEV root 2026110701 (+ shard). Output: deposits/dt8/dev_n{n}_s{shard}.json. Not a CONFIRM run.
Usage: python dt8_pinned_dev.py n R shard nshards
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pinned_dip import Y_joined, Z_dec  # noqa: E402

US = [0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.93, 0.95, 0.97, 0.99, 1.0]


def main():
    n, R, shard, ns = (int(a) for a in sys.argv[1:5])
    out = os.path.join(os.path.dirname(__file__), "..", "deposits", "dt8")
    os.makedirs(out, exist_ok=True)
    ss = np.random.SeedSequence([2026110701, n, shard])
    rng = np.random.default_rng(ss)
    Ys, Zs = [], []
    t0 = time.time()
    for _ in range(shard, R, ns):
        x = np.sort(rng.uniform(size=n))
        Ys.append([np.sqrt(n) * Y_joined(x, u) for u in US])
        Zs.append([np.sqrt(n) * Z_dec(x, u) for u in US])
    res = {"n": n, "R_shard": len(Ys), "shard": shard, "nshards": ns, "us": US,
           "Y": Ys, "Z": Zs, "seconds": time.time() - t0, "root": 2026110701}
    with open(os.path.join(out, f"dev_n{n}_s{shard}.json"), "w") as f:
        json.dump(res, f)
    print(f"n={n} shard={shard} reps={len(Ys)} {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
