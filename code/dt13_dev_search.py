"""Step 14 (DT13) DEV search (root 2026110743): per-n worst two-step law W(w, eps) for the uniform-calibrated dip.

Family (CONVENTIONS CV-9; DT2 K1/K7: one flat shoulder + a top step is the strongest structure found):
  W(w, eps): density (1-w)/L on (0, L), w/eps on (L, 1), L = 1 - eps.  (Unimodal for every w, eps.)
Grid: w in WS, eps in ES. Per n: DEV uniform reference R_REF (spawn (0, n)), each law N draws (spawn (1, n, iw, ie)).
Selection is on DEV noise (winner's curse, N-2): the chosen laws are re-estimated on a FRESH CONFIRM root later.
Usage: python dt13_dev_search.py  -> deposits/dt13/dev_search.json
"""
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dt0c_screen import dips_of  # noqa: E402

ROOT = 2026110743
NS = (5, 6, 7, 8, 10, 12, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 300, 500, 750, 1000)
WS = (0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7)
ES = (1e-6, 1e-4, 1e-3, 0.005, 0.02)
ALPHAS = (0.01, 0.05, 0.10)
R_REF, N = 200_000, 20_000


def draw_w(rng, size, w, eps):
    L = 1.0 - eps
    u = rng.random(size)
    top = rng.random(size) < w
    return np.where(top, L + eps * u, L * u)


def batches(rng, m, n, gen):
    b = max(200, min(20000, 2_000_000 // n))
    out = np.empty(m)
    for s in range(0, m, b):
        k = min(b, m - s)
        out[s:s + k] = dips_of(gen(rng, (k, n)))
    return out


def job(n):
    t0 = time.time()
    rng = np.random.default_rng(np.random.SeedSequence([ROOT, 0, n]))
    ref = np.sort(batches(rng, R_REF, n, lambda r, s: r.random(s)))
    q = {a: float(ref[int(np.ceil((1 - a) * R_REF)) - 1]) for a in ALPHAS}
    rows = []
    for iw, w in enumerate(WS):
        for ie, e in enumerate(ES):
            rg = np.random.default_rng(np.random.SeedSequence([ROOT, 1, n, iw, ie]))
            d = batches(rg, N, n, lambda r, s, w=w, e=e: draw_w(r, s, w, e))
            rows.append({"w": w, "eps": e, **{str(a): float(np.mean(d >= q[a])) for a in ALPHAS}})
    best = {}
    for a in ALPHAS:
        b = max(rows, key=lambda r, a=a: r[str(a)])
        best[str(a)] = b
    return n, {"q": q, "rows": rows, "best": best, "seconds": time.time() - t0}


def main():
    res = {"root": ROOT, "R_ref": R_REF, "N": N, "WS": WS, "ES": ES, "by_n": {}}
    with ProcessPoolExecutor(max_workers=11) as ex:
        for n, r in ex.map(job, sorted(NS, key=lambda v: -v)):
            res["by_n"][str(n)] = r
            print(f"n={n:5d} {r['seconds']:.0f}s " + " | ".join(
                f"a={a}: w={r['best'][str(a)]['w']} e={r['best'][str(a)]['eps']:.0e} {100 * r['best'][str(a)][str(a)]:.2f}%"
                for a in ALPHAS), flush=True)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "deposits", "dt13")
    try:
        os.makedirs(out, exist_ok=True)
        with open(os.path.join(out, "dev_search.json"), "w") as fh:
            json.dump(res, fh, indent=1)
    except OSError as err:
        sys.exit(f"cannot write: {err}")


if __name__ == "__main__":
    main()
