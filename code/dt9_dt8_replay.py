"""DT9 Part R: replay of the DT8 Part B CONFIRM samples (root 2026110711) with certified numerics.

DT8 used LP values without a repair term (its 1e-9 margin was not proved to cover shape violations).
This script regenerates EXACTLY the DT8 samples (same SeedSequence([root, n, shard]) stream, same order:
for _ in range(shard, R, ns): x = sort(uniform(n))) and recomputes, per block, the conservative indicators
of dt9_block_cert (certified primal upper bounds; failure -> exceedance) at every DT8 grid cutoff.
Cutoffs are evaluated in increasing order; a higher one only after an exceedance at the lower one.
No new random draws: this can only confirm or withdraw DT8 claims, never add new ones.
Usage: python dt9_dt8_replay.py n R shard nshards   -> deposits/dt9/replay_n<n>_s<shard>.json
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dt9_block_cert as C  # noqa: E402
from pinned_dip import V_left  # noqa: E402
from sup_pinned import pieces  # noqa: E402

ROOT = 2026110711
CUTS = (0.7315, 0.735, 0.74, 0.75, 0.86, 0.87, 0.88, 0.90)   # increasing; both DT8 grids


def main():
    n, R, shard, ns = (int(a) for a in sys.argv[1:5])
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "deposits", "dt9")
    path = os.path.join(out, f"replay_n{n}_s{shard}.json")
    if os.path.exists(path):
        return
    rng = np.random.default_rng(np.random.SeedSequence([ROOT, n, shard]))
    rt = np.sqrt(n)
    t0 = time.time()
    EX, FL, VV = [], [], []
    for _ in range(shard, R, ns):
        x = np.sort(rng.uniform(size=n))
        VV.append(float(rt * V_left(x, 1.0)))   # identity check: DT8 stored V for EVERY sample
        row, fl = [], []
        for lo_u, hi_u in C.BLOCKS:
            P = pieces(x, lo_u, hi_u)
            ex, anyf = [], False
            for c in CUTS:
                if ex and not ex[-1]:
                    ex.append(False)
                    continue
                e, f = C.exceeds(x, P, c / rt)
                anyf = anyf or f
                ex.append(bool(e))
            row.append(ex)
            fl.append(anyf)
        EX.append(row)
        FL.append(fl)
    res = {"n": n, "root": ROOT, "shard": shard, "nshards": ns, "R": R, "cuts": CUTS, "ex": EX, "failed": FL, "V": VV,
           "stats": dict(C.STATS), "seconds": time.time() - t0}
    try:
        with open(path + ".tmp", "w") as fh:
            json.dump(res, fh)
        os.replace(path + ".tmp", path)
    except OSError as err:
        sys.exit(f"cannot write {path}: {err}")
    print(f"replay n={n} shard={shard} reps={len(EX)} {time.time() - t0:.1f}s lps={C.STATS['lps']} "
          f"fb={C.STATS['fallbacks']} rep={C.STATS['max_repair']:.1e}", flush=True)


if __name__ == "__main__":
    main()
