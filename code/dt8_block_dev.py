"""DT8 DEV pilot (Step 10; DT7-UNIF Thm U1): block envelopes Sbar_I for uniform samples.

For each uniform sample of size n and each block I of pinned modes in [1/2, 1], computes
sqrt(n) * Sbar_I (branch and bound; exact when >= FLOOR, else an upper bound < FLOOR), plus the
edge statistic sqrt(n) V (u = 1). Used ONLY to choose the block partition and the cutoff that the
frozen CONFIRM run will test. DEV root 2026110703 (+ n, shard). Not a CONFIRM run.
Usage: python dt8_block_dev.py n R shard nshards [root] [tag]
CONFIRM use (DT8-PREREG part B): tag = cert; a block whose LP fails after all fallbacks is recorded as +inf
(counted as an exceedance at every cutoff; conservative policy) and listed in "failed".
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pinned_dip import V_left  # noqa: E402
from sup_pinned import FALLBACKS, MAXVIOL, block_stat_bb  # noqa: E402

BLOCKS = [(0.5, 0.9), (0.9, 0.97), (0.97, 0.995), (0.995, 1.0)]
FLOOR = 0.65          # sqrt(n) scale


def main():
    n, R, shard, ns = (int(a) for a in sys.argv[1:5])
    root = int(sys.argv[5]) if len(sys.argv) > 5 else 2026110703
    tag = sys.argv[6] if len(sys.argv) > 6 else "dev"
    failed = []
    out = os.path.join(os.path.dirname(__file__), "..", "deposits", "dt8blk")
    os.makedirs(out, exist_ok=True)
    rng = np.random.default_rng(np.random.SeedSequence([root, n, shard]))
    rt = np.sqrt(n)
    S, X, V = [], [], []
    t0 = time.time()
    for _ in range(shard, R, ns):
        x = np.sort(rng.uniform(size=n))
        row, ex = [], []
        for jb, (lo_u, hi_u) in enumerate(BLOCKS):
            try:
                v, e = block_stat_bb(x, lo_u, hi_u, FLOOR / rt)
            except RuntimeError:
                v, e = float("inf"), False
                failed.append([len(S), jb])
            row.append(rt * v)
            ex.append(bool(e))
        S.append(row)
        X.append(ex)
        V.append(rt * V_left(x, 1.0))
    res = {"n": n, "root": root, "shard": shard, "nshards": ns, "blocks": BLOCKS, "floor": FLOOR,
           "S": S, "exact": X, "V": V, "maxviol": MAXVIOL[0], "fallbacks": FALLBACKS[0], "failed": failed,
           "seconds": time.time() - t0}
    with open(os.path.join(out, f"{tag}_n{n}_r{root}_s{shard}.json"), "w") as f:
        json.dump(res, f)
    print(f"n={n} shard={shard} reps={len(S)} {time.time() - t0:.1f}s maxviol={MAXVIOL[0]:.1e} "
          f"fallbacks={FALLBACKS[0]}")


if __name__ == "__main__":
    main()
