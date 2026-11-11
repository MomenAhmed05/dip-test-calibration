"""Summarise DT8 DEV pinned-envelope shards: P(sqrt(n) Y(u) > t), P(sqrt(n) Z(u) > t), with 95% CP intervals."""
import glob
import json
import os
import sys

import numpy as np
from scipy.stats import beta

TS = [0.5452, 0.73139, 0.85987]


def cp(k, N, a=0.05):
    lo = 0.0 if k == 0 else beta.ppf(a / 2, k, N - k + 1)
    hi = 1.0 if k == N else beta.ppf(1 - a / 2, k + 1, N - k)
    return lo, hi


def main():
    d = os.path.join(os.path.dirname(__file__), "..", "deposits", "dt8")
    for n in sys.argv[1:]:
        files = sorted(glob.glob(os.path.join(d, f"dev_n{n}_s*.json")))
        if not files:
            continue
        Y, Z, us = [], [], None
        for f in files:
            r = json.load(open(f))
            Y += r["Y"]; Z += r["Z"]; us = r["us"]
        Y, Z = np.array(Y), np.array(Z)
        N = len(Y)
        print(f"n={n} N={N} shards={len(files)}")
        for t in TS:
            print(f"  t={t}")
            for j, u in enumerate(us):
                ky, kz = int((Y[:, j] > t).sum()), int((Z[:, j] > t).sum())
                ly, hy = cp(ky, N)
                print(f"    u={u:5.2f}  Y {ky / N:7.4f} [{ly:.4f},{hy:.4f}]   Z {kz / N:7.4f}")
        # pathwise: Y(u) <= Y(1)?  fraction of paths where interior Y exceeds edge
        ex = (Y[:, :-1] > Y[:, -1:] + 1e-12).mean(axis=0)
        print("  frac paths Y(u) > Y(1):", np.round(ex, 3).tolist())


if __name__ == "__main__":
    main()
