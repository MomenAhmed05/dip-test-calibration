"""DT10 frozen analysis (DT10-PREREG): empirical-Bernstein upper bounds for E[V_i] per eps-interval and level.

Maurer & Pontil (2009), Thm 4: for i.i.d. V in [0, 1], n >= 2, with probability >= 1 - d,
    E V <= mean + sqrt(2 s^2 ln(2/d) / n) + 7 ln(2/d) / (3 (n - 1)),   s^2 = unbiased sample variance.
d = DELTA / (number of claims) (Bonferroni over all intervals x levels; shared paths are allowed).
Interval i at level alpha is CERTIFIED iff its bound is < alpha. Reported range: the union of certified
intervals, and the largest contiguous certified block that ends at eps = 1/sqrt(2).
Usage: python dt10_cert_analysis.py root tag m  -> deposits/dt10/cert_analysis_<tag>_r<root>.json
"""
import glob
import json
import math
import os
import sys

import numpy as np

DELTA = 0.01
NSH = 11


def main():
    root, tag, m = int(sys.argv[1]), sys.argv[2], int(sys.argv[3])
    base = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "deposits", "dt10")
    fs = sorted(glob.glob(os.path.join(base, f"cert_{tag}_m{m}_r{root}_s*.json")))
    if len(fs) != NSH:
        sys.exit(f"expected {NSH} shard files, found {len(fs)}")
    shards = []
    for f in fs:
        try:
            with open(f) as fh:
                shards.append(json.load(fh))
        except (OSError, ValueError) as err:
            sys.exit(f"cannot read {f}: {err}")
    R, edges, levels = shards[0]["R"], list(shards[0]["edges"]), list(shards[0]["levels"])
    n = sum(int(d["reps"]) for d in shards)
    S1 = np.sum([np.array(d["S1"], float) for d in shards], axis=0)
    S2 = np.sum([np.array(d["S2"], float) for d in shards], axis=0)
    if any(d["edges"] != edges or d["levels"] != levels or d["R"] != R for d in shards):
        sys.exit("shard metadata differ")
    if n != R:
        sys.exit(f"replication count {n} != planned {R}")
    mean = S1 / n
    var = np.maximum((S2 - n * mean * mean) / (n - 1), 0.0)
    claims = mean.size
    L = math.log(2.0 * claims / DELTA)
    ub = mean + np.sqrt(2.0 * var * L / n) + 7.0 * L / (3.0 * (n - 1))
    rows = []
    for li, (alpha, t_lo) in enumerate(levels):
        cert = [bool(ub[li, i] + 1e-12 < alpha) for i in range(len(edges) - 1)]   # pad covers float rounding
        k = len(cert)
        while k > 0 and cert[k - 1]:
            k -= 1
        rows.append({"alpha": alpha, "t_lo": t_lo, "mean": mean[li].tolist(), "sd": np.sqrt(var[li]).tolist(),
                     "upper": ub[li].tolist(), "certified": cert,
                     "contiguous_eps_from": edges[k] if k < len(cert) else None,
                     "contiguous_u_upto": (1.0 - edges[k] ** 2) if k < len(cert) else None,
                     "n_certified": int(sum(cert))})
        print(f"alpha={alpha} t_lo={t_lo} n={n} claims={claims} certified {sum(cert)}/{len(cert)}; contiguous "
              f"eps in [{rows[-1]['contiguous_eps_from']}, 1/sqrt2]; min margin {float(np.min(alpha - ub[li])):+.5f}")
    out = os.path.join(base, f"cert_analysis_{tag}_r{root}.json")
    try:
        with open(out, "w") as fh:
            json.dump({"root": root, "tag": tag, "m": m, "n": n, "delta": DELTA, "claims": claims, "L": L,
                       "edges": edges, "rows": rows}, fh, indent=1)
    except OSError as err:
        sys.exit(f"cannot write {out}: {err}")


if __name__ == "__main__":
    main()
