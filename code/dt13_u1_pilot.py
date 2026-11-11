"""Step 14 (DT13) DEV pilot (root 2026110742): is the Corollary U1' block bound informative at the UNIFORM cutoff?

For a uniform sample of size n and a block I of pinned modes, the conservative indicator Z_I(t) = 1 unless a
certified upper bound shows Sbar_I < t (dt9_block_cert.exceeds with a strict threshold t - 1e-12). We estimate
P(Z_I(t) = 1) for t = q_n(1 - alpha) (raw dip scale), alpha in {.01, .05, .10}, for two block partitions:
  B4 = DT8/DT9 blocks [.5,.9], [.9,.97], [.97,.995], [.995,1]
  B8 = a finer partition [.5,.7], [.7,.8], [.8,.9], [.9,.95], [.95,.97], [.97,.985], [.985,.995], [.995,1]
Cutoffs q_n: n = 4 exact (dt13_n4_exact); other n from a quick DEV uniform reference (R_ref per n, same root).
Usage: python dt13_u1_pilot.py n N [Rref] ; writes deposits/dt13/pilot_u1_n{n}.json
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dt0c_screen import dips_of  # noqa: E402
from dt9_block_cert import exceeds  # noqa: E402
from sup_pinned import pieces  # noqa: E402

ROOT = 2026110742
B4 = [(0.5, 0.9), (0.9, 0.97), (0.97, 0.995), (0.995, 1.0)]
B8 = [(0.5, 0.7), (0.7, 0.8), (0.8, 0.9), (0.9, 0.95), (0.95, 0.97), (0.97, 0.985), (0.985, 0.995), (0.995, 1.0)]
ALPHAS = (0.01, 0.05, 0.10)
Q4 = {0.01: None, 0.05: None, 0.10: 3 / 16}   # n = 4: q = r/(4(r+1)); r_.10 = 3 -> 3/16


def q4(a):
    r = (np.sqrt(1 + 8 / a) - 3) / 2
    return r / (4 * (r + 1))


def main():
    n, N = int(sys.argv[1]), int(sys.argv[2])
    Rref = int(sys.argv[3]) if len(sys.argv) > 3 else 200_000
    rng = np.random.default_rng(np.random.SeedSequence([ROOT, n]))
    if n == 4:
        qs = {a: q4(a) for a in ALPHAS}
    else:
        ref = np.sort(np.concatenate([dips_of(rng.random((min(20000, Rref), n))) for _ in range(max(1, Rref // 20000))]))
        qs = {a: float(ref[int(np.ceil((1 - a) * len(ref))) - 1]) for a in ALPHAS}
    t0 = time.time()
    cnt = {nm: np.zeros((len(B), len(ALPHAS)), int) for nm, B in (("B4", B4), ("B8", B8))}
    for _ in range(N):
        x = np.sort(rng.random(n))
        for nm, B in (("B4", B4), ("B8", B8)):
            for j, (lo, hi) in enumerate(B):
                P = pieces(x, lo, hi)
                prev = True
                for ia, a in sorted(enumerate(ALPHAS), key=lambda z: -z[1]):   # largest alpha = smallest cut first
                    if not prev:          # nested: higher cutoff cannot be exceeded if a lower one is not
                        continue
                    e, _f = exceeds(x, P, qs[a] - 1e-12)
                    cnt[nm][j, ia] += int(e)
                    prev = e
    res = {"n": n, "N": N, "q": qs, "sqrt_n_q": {a: float(np.sqrt(n) * q) for a, q in qs.items()},
           "seconds": time.time() - t0}
    for nm, B in (("B4", B4), ("B8", B8)):
        res[nm] = {"blocks": B, "counts": cnt[nm].tolist(),
                   "rates": (cnt[nm] / N).tolist(), "max_rate": (cnt[nm] / N).max(axis=0).tolist(),
                   "argmax_block": cnt[nm].argmax(axis=0).tolist()}
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "deposits", "dt13")
    try:
        os.makedirs(out, exist_ok=True)
        with open(os.path.join(out, f"pilot_u1_n{n}.json"), "w") as fh:
            json.dump(res, fh, indent=1)
    except OSError as err:
        sys.exit(f"cannot write pilot output: {err}")
    print(f"n={n} N={N} {res['seconds']:.0f}s sqrt(n)q={[round(v, 4) for v in res['sqrt_n_q'].values()]} "
          f"B4 max={[round(v, 4) for v in res['B4']['max_rate']]} at {res['B4']['argmax_block']} "
          f"B8 max={[round(v, 4) for v in res['B8']['max_rate']]} at {res['B8']['argmax_block']}", flush=True)


if __name__ == "__main__":
    main()
